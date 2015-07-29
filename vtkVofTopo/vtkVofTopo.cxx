#include "vtkVofTopo.h"
#include "vofTopology.h"

#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkMPIController.h"
#include "vtkMPICommunicator.h"
#include "vtkPolyData.h"
#include <iostream>
#include <vector>
#include <cmath>

vtkStandardNewMacro(vtkVofTopo);

// taken from vtkParticleTracerBase.cxx
namespace
{
  int findClosestTimeStep(double requestedTimeValue,
			  const std::vector<double>& timeSteps)
  {
    int ts = 0;
    double mindist = std::abs(timeSteps[0] - requestedTimeValue);

    for (int i = 0; i < timeSteps.size(); i++) {

      double tsv = timeSteps[i];
      double dist = std::abs(tsv - requestedTimeValue);
      if (dist < mindist) {
	mindist = dist;
	ts = i;
      }
    }
    return ts;
  }
};

//----------------------------------------------------------------------------
vtkVofTopo::vtkVofTopo() :
  FirstIteration(true),
  CurrentTimeStep(0),
  LastComputedTimeStep(-1),
  UseCache(false),
  IterType(IterateOverTarget),
  Seeds(0)
{
  this->SetNumberOfInputPorts(2);
  this->Controller = vtkMPIController::New();
}

//----------------------------------------------------------------------------
vtkVofTopo::~vtkVofTopo()
{
  if (Seeds != 0) {
    Seeds->Delete();
  }
  this->Controller->Delete();
}

//----------------------------------------------------------------------------
int vtkVofTopo::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  }
  if (port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkVofTopo::RequestInformation(vtkInformation *vtkNotUsed(request),
				   vtkInformationVector **inputVector,
				   vtkInformationVector *outputVector)
{
  vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);

  if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {

    unsigned int numberOfInputTimeSteps =
      inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

    this->InputTimeValues.resize(numberOfInputTimeSteps);
    inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
		&this->InputTimeValues[0]);

    if (numberOfInputTimeSteps == 1) {
      vtkWarningMacro(<<"Not enough input time steps for topology computation");
    }

    // check if user input is within time step range
    if (InitTimeStep < 0) {
      InitTimeStep = 0;
    }
    if (InitTimeStep > InputTimeValues.size()-1) {
      InitTimeStep = InputTimeValues.size()-1;
    }
    if (TargetTimeStep > InputTimeValues.size()-1) {
      TargetTimeStep = InputTimeValues.size()-1;
    }
    if (TargetTimeStep < InitTimeStep) {
      return 0;
    }

    // multiprocess
    if (Controller->GetCommunicator() != 0) {
      // find neighbor processes and global domain bounds
      GetGlobalContext(inInfo);
    }
  }
  else {
    vtkErrorMacro(<<"Input information has no TIME_STEPS set");
    return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkVofTopo::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
				    vtkInformationVector **inputVector,
				    vtkInformationVector *outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  // inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  // outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);

  if(FirstIteration) {

    if (IterType == IterateOverTarget) {

      double targetTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      if(targetTime > InputTimeValues.back()) {
      	targetTime = InputTimeValues.back();
      }
      TargetTimeStep = findClosestTimeStep(targetTime, InputTimeValues);
      if (TargetTimeStep < 0) {
	TargetTimeStep = 0;
      }
      if (TargetTimeStep > InputTimeValues.size()-1) {
	TargetTimeStep = InputTimeValues.size()-1;
      }

      if (LastComputedTimeStep > -1 &&
      	  LastComputedTimeStep < TargetTimeStep) {
      	UseCache = true;
      }
      else {
      	UseCache = false;
      	LastComputedTimeStep = -1;
      }
      
      if (UseCache) {
        CurrentTimeStep = LastComputedTimeStep + 1;
      }
      else {
	CurrentTimeStep = InitTimeStep;
      }
    }
  }
  if (CurrentTimeStep <= TargetTimeStep) {
    std::cout << "CurrentTimeStep = " << CurrentTimeStep << std::endl;
    int numInputs = this->GetNumberOfInputPorts();
    for (int i = 0; i < numInputs; i++) {
      vtkInformation *inInfo = inputVector[i]->GetInformationObject(0);

      if (CurrentTimeStep < static_cast<int>(InputTimeValues.size())) {
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(),
		    InputTimeValues[CurrentTimeStep]);
      }
    }
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkVofTopo::RequestData(vtkInformation *request,
			    vtkInformationVector **inputVector,
			    vtkInformationVector *outputVector)
{
  vtkInformation *inInfoVelocity = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfoVof = inputVector[1]->GetInformationObject(0);

  vtkRectilinearGrid *inputVelocity = vtkRectilinearGrid::
    SafeDownCast(inInfoVelocity->Get(vtkDataObject::DATA_OBJECT()));
  vtkRectilinearGrid *inputVof = vtkRectilinearGrid::
    SafeDownCast(inInfoVof->Get(vtkDataObject::DATA_OBJECT()));

  if (IterType == IterateOverTarget) {

    // Stage I ---------------------------------------------------------------
    if (FirstIteration) {
      if (!UseCache) {
	GenerateSeeds(inputVof);
	InitParticles();
      }
    }

    // Stage II --------------------------------------------------------------
    if(CurrentTimeStep < TargetTimeStep) {
      AdvectParticles(inputVof, inputVelocity);
      LastComputedTimeStep = CurrentTimeStep;
    }

    bool finishedAdvection = CurrentTimeStep >= TargetTimeStep;
    if (finishedAdvection) {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      FirstIteration = true;

      vtkSmartPointer<vtkRectilinearGrid> components = vtkSmartPointer<vtkRectilinearGrid>::New();
      // Stage III -----------------------------------------------------------
      ExtractComponents(inputVof, components);

      vtkInformation *outInfo = outputVector->GetInformationObject(0);
      vtkRectilinearGrid *output =
      	vtkRectilinearGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

      output->SetExtent(components->GetExtent());
      output->SetXCoordinates(components->GetXCoordinates());
      output->SetYCoordinates(components->GetYCoordinates());
      output->SetZCoordinates(components->GetZCoordinates());
      output->GetCellData()->AddArray(components->GetCellData()->GetArray("Labels"));
      output->GetCellData()->SetActiveScalars("Labels");

      // vtkInformation *outInfo = outputVector->GetInformationObject(0);
      // vtkPolyData *output =
      // 	vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
      // GenerateOutputGeometry(output);
    }
    else {
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
      FirstIteration = false;
      CurrentTimeStep++;
    }
  }
  else { // IterType == IterateOverInit

  }

  return 1;
}

//----------------------------------------------------------------------------
void vtkVofTopo::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkVofTopo::GenerateSeeds(vtkRectilinearGrid *vof)
{
  vtkSmartPointer<vtkPoints> seedPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> seedConnectivity = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkShortArray> seedCoords = vtkSmartPointer<vtkShortArray>::New();
  generateSeedPoints(vof, Refinement, seedPoints, seedConnectivity, seedCoords);
  if (Seeds != 0) {
    Seeds->Delete();
  }
  Seeds = vtkPolyData::New();
  Seeds->SetPoints(seedPoints);
  Seeds->GetPointData()->AddArray(seedConnectivity);
  Seeds->GetPointData()->AddArray(seedCoords);
}

//----------------------------------------------------------------------------
void vtkVofTopo::InitParticles()
{
  vtkPoints *seedPoints = Seeds->GetPoints();
  Particles.clear();
  for (int i = 0; i < seedPoints->GetNumberOfPoints(); ++i) {
    double p[3];
    seedPoints->GetPoint(i, p);
    Particles.push_back(make_float4(p[0], p[1], p[2], 1.0f));
  }
}

//----------------------------------------------------------------------------
void vtkVofTopo::GenerateOutputGeometry(vtkPolyData *output)
{
  vtkPoints *advectedParticles = vtkPoints::New();
  for (int i = 0; i < Particles.size(); ++i) {
    double p[3] = {Particles[i].x, Particles[i].y, Particles[i].z};
    advectedParticles->InsertNextPoint(p);
  }
  output->SetPoints(advectedParticles);
}

//----------------------------------------------------------------------------
void vtkVofTopo::GetGlobalContext(vtkInformation *inInfo)
{
  int numProcesses = Controller->GetNumberOfProcesses();
  std::vector<vtkIdType> RecvLengths(numProcesses);
  std::vector<vtkIdType> RecvOffsets(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    RecvLengths[i] = 6;
    RecvOffsets[i] = i*6;
  }

  vtkRectilinearGrid *inputVof = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
      
  int LocalExtent[6];
  inputVof->GetExtent(LocalExtent);
  std::vector<int> AllExtents(6*numProcesses);
  Controller->AllGatherV(&LocalExtent[0], &AllExtents[0], 6, &RecvLengths[0], &RecvOffsets[0]);

  int GlobalExtents[6];
  findGlobalExtents(AllExtents, GlobalExtents);

  NeighborProcesses.clear();
  NeighborProcesses.resize(6);
  findNeighbors(LocalExtent, GlobalExtents, AllExtents, NeighborProcesses);

  NumNeighbors = 0;
  for (int i = 0; i < NeighborProcesses.size(); ++i) {
    NumNeighbors += NeighborProcesses[i].size();
  }

  // find domain bounds
  inputVof->GetBounds(&LocalBounds[0]);
  std::vector<double> AllBounds(6*numProcesses);
  Controller->AllGatherV(&LocalBounds[0], &AllBounds[0], 6, &RecvLengths[0], &RecvOffsets[0]);
  findGlobalBounds(AllBounds, GlobalBounds);
}

//----------------------------------------------------------------------------
void vtkVofTopo::AdvectParticles(vtkRectilinearGrid *vof,
				 vtkRectilinearGrid *velocity)
{
  float dt = InputTimeValues[CurrentTimeStep+1] - InputTimeValues[CurrentTimeStep];
  advectParticles(vof, velocity, Particles, dt);
  if (Controller->GetCommunicator() != 0) {
    ExchangeParticles();
  }
}

//----------------------------------------------------------------------------
void vtkVofTopo::ExchangeParticles()
{
  int numProcesses = Controller->GetNumberOfProcesses();
  int processId = Controller->GetLocalProcessId();

  // one vector for each side of the process
  std::vector<std::vector<float4> > particlesToSend(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    particlesToSend[i].resize(0);
  }

  std::vector<float4>::iterator it;
  std::vector<float4> particlesToKeep;

  for (it = Particles.begin(); it != Particles.end(); ++it) {

    int bound = outOfBounds(*it, LocalBounds, GlobalBounds);
    if (bound > -1) {
      for (int j = 0; j < NeighborProcesses[bound].size(); ++j) {

  	int neighborId = NeighborProcesses[bound][j];	
  	particlesToSend[neighborId].push_back(*it);
      }
    }
    else {
      particlesToKeep.push_back(*it);
    }
  }
  Particles = particlesToKeep;

  std::vector<float4> particlesToRecv;
  sendData(particlesToSend, particlesToRecv, numProcesses, Controller);

  // insert the paricles that are within the domain
  for (int i = 0; i < particlesToRecv.size(); ++i) {
    int within = withinBounds(particlesToRecv[i], LocalBounds);
    if (within) {
      Particles.push_back(particlesToRecv[i]);
    }
  }
}

//----------------------------------------------------------------------------
void vtkVofTopo::ExtractComponents(vtkRectilinearGrid *vof,
				   vtkRectilinearGrid *components)
{
  int nodeRes[3];
  vof->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

  vtkDataArray *data = vof->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

  vtkFloatArray *labels = vtkFloatArray::New();
  labels->SetName("Labels");
  labels->SetNumberOfComponents(1);
  labels->SetNumberOfTuples(vof->GetNumberOfCells());
  for (int i = 0; i < labels->GetNumberOfTuples(); ++i) {
    labels->SetValue(i, -1.0f);
  }

  g_labelId = 0;
  // determine if data is float or double
  if (data->IsA("vtkFloatArray")) {
    extractComponents(vtkFloatArray::SafeDownCast(data)->GetPointer(0),
  		      cellRes, labels->GetPointer(0));
  }
  else if (data->IsA("vtkDoubleArray")) {
    extractComponents(vtkDoubleArray::SafeDownCast(data)->GetPointer(0),
  		      cellRes, labels->GetPointer(0));
  }
   
  //--------------------------------------------------------------------------
  // send number of labels to other processes
  if (Controller->GetCommunicator() != 0) {

    int numProcesses = this->Controller->GetNumberOfProcesses();
    int processId = Controller->GetLocalProcessId();

    // -----------------------------------------------------------------------
    // gather number of labels from other processes
    std::vector<vtkIdType> recvLengths(numProcesses);
    std::vector<vtkIdType> recvOffsets(numProcesses);
    for (int i = 0; i < numProcesses; ++i) {
      recvLengths[i] = 1;
      recvOffsets[i] = i;
    }

    int numMyLabels = g_labelId;
    std::vector<int> allNumLabels(numProcesses);
    Controller->AllGatherV(&numMyLabels, &allNumLabels[0], 1, &recvLengths[0], &recvOffsets[0]);
    std::vector<int> labelOffsets(numProcesses);
    labelOffsets[0] = 0;
    for (int i = 1; i < numProcesses; ++i) {
      labelOffsets[i] = labelOffsets[i-1] + allNumLabels[i-1];
    }
    int numAllLabels = labelOffsets.back() + allNumLabels.back();

    // ------------------------
    for (int i = 0; i < numProcesses; ++i) {
      recvLengths[i] = 6;
      recvOffsets[i] = i*6;
    }

    // int MyExtent[6];
    // input->GetExtent(MyExtent);
    // std::vector<int> AllExtents(6*numProcesses);
    // Controller->AllGatherV(&MyExtent[0], &AllExtents[0], 6, &recvLengths[0], &recvOffsets[0]);

    // int GlobalExtents[6];
    // findGlobalExtents(AllExtents, GlobalExtents);
    // std::vector<std::vector<int> > neighborProcesses;
    // neighborProcesses.clear();
    // neighborProcesses.resize(6);
    // findNeighbors(MyExtent, GlobalExtents, AllExtents, neighborProcesses);

    // int numNeighbors = 0;
    // for (int i = 0; i < neighborProcesses.size(); ++i) {
    //   numNeighbors += neighborProcesses[i].size();
    // }

    // -----------------------------------------------------------------------
    // prepare labelled cells to send to neighbors
    int myExtent[6];
    vof->GetExtent(myExtent);

    int slabs[6][6] = {{0,0,                       0,cellRes[1]-1,            0,cellRes[2]-1},
		       {cellRes[0]-1,cellRes[0]-1, 0,cellRes[1]-1,            0,cellRes[2]-1},
		       {0,cellRes[0]-1,            0,0,                       0,cellRes[2]-1},
		       {0,cellRes[0]-1,            cellRes[1]-1,cellRes[1]-1, 0,cellRes[2]-1},
		       {0,cellRes[0]-1,            0,cellRes[1]-1,            0,0},
		       {0,cellRes[0]-1,            0,cellRes[1]-1,            cellRes[2]-1,cellRes[2]-1}};

    std::vector<std::vector<float4> > labelsToSend(6);
    for (int p = 0; p < NeighborProcesses.size(); ++p) {
      if (NeighborProcesses[p].size() > 0) {

    	for (int k = slabs[p][4]; k <= slabs[p][5]; ++k) {
    	  for (int j = slabs[p][2]; j <= slabs[p][3]; ++j) {
    	    for (int i = slabs[p][0]; i <= slabs[p][1]; ++i) {

    	      int idx = i + j*cellRes[0] + k*cellRes[0]*cellRes[1];
    	      float label = labels->GetValue(idx);
    	      if (label > -1) {
    		labelsToSend[p].push_back(make_float4(i+myExtent[0], j+myExtent[2],
    						      k+myExtent[4], label));
    	      }
    	    }
    	  }
    	}
      }
    }

    // -----------------------------------------------------------------------
    // send header to neighbors with the number of labels to be send
    int numLabelsToSend[6];
    for (int i = 0; i < 6; ++i) {
      numLabelsToSend[i] = labelsToSend[i].size();

      for (int j = 0; j < NeighborProcesses[i].size(); ++j) {

    	const int SEND_LABELS_TAG = 100+processId;
    	Controller->Send(&numLabelsToSend[i], 1,
    			 NeighborProcesses[i][j], SEND_LABELS_TAG);
      }
    }

    // -----------------------------------------------------------------------
    // receive header
    std::vector<int> numLabelsToRecv(NumNeighbors);
    int nidx = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < NeighborProcesses[i].size(); ++j) {
    	numLabelsToRecv[nidx] = 0;
    	const int RECV_LABELS_TAG = 100+NeighborProcesses[i][j];
    	Controller->Receive(&numLabelsToRecv[nidx], 1,
    			    NeighborProcesses[i][j], RECV_LABELS_TAG);
    	++nidx;
      }
    }

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------

    // send the labels to each side
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < NeighborProcesses[i].size(); ++j) {
    	const int SEND_LABEL_DATA_TAG = 100+processId;

    	vtkMPICommunicator::Request req;
    	Controller->NoBlockSend((char*)&(labelsToSend[i][0]),
    				labelsToSend[i].size()*sizeof(float4),
    				NeighborProcesses[i][j], SEND_LABEL_DATA_TAG, req);
      }
    }
    // -----------------------------------------------------------------------
    // allocate buffers to receive labels from each neighbor
    std::vector<std::vector<float4> > labelsToRecv(NumNeighbors);
    for (int i = 0; i < NumNeighbors; ++i) {
      labelsToRecv[i].resize(numLabelsToRecv[i]);
    }
    // -----------------------------------------------------------------------
    // receive labels from each neighbor
    vtkMPICommunicator::Request *reqs = new vtkMPICommunicator::Request[NumNeighbors];
    nidx = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < NeighborProcesses[i].size(); ++j) {
    	const int RECV_LABEL_DATA_TAG = 100+NeighborProcesses[i][j];

    	Controller->NoBlockReceive((char*)&(labelsToRecv[nidx][0]),
    				   labelsToRecv[nidx].size()*sizeof(float4),
    				   NeighborProcesses[i][j], RECV_LABEL_DATA_TAG, reqs[nidx]);
    	++nidx;
      }
    }
    Controller->WaitAll(NumNeighbors, reqs);

    std::vector<int> allLabels(numAllLabels);
    for (int i = 0; i < allLabels.size(); ++i) {
      allLabels[i] = i;
    }
    
    // // -----------------------------------------------------------------------
    // // identify equivalent labels from neighbor processes
    // nidx = 0;
    // for (int i = 0; i < 6; ++i) {
    //   for (int j = 0; j < NeighborProcesses[i].size(); ++j) {
    // 	for (int s = 0; s < labelsToRecv[nidx].size(); ++s) {

    // 	  int x = labelsToRecv[nidx][s].x;
    // 	  int y = labelsToRecv[nidx][s].y;
    // 	  int z = labelsToRecv[nidx][s].z;
    // 	  int neighborLabel = labelsToRecv[nidx][s].w + labelOffsets[NeighborProcesses[i][j]];
	  
    // 	  int mx, my, mz;

    // 	  if (!adjointPoint(x, y, z, myExtent)) {
    // 	    continue;
    // 	  }
    // 	  mx = x - myExtent[0];
    // 	  my = y - myExtent[2];
    // 	  mz = z - myExtent[4];

    // 	  int idx = mx + my*res[0] + mz*res[0]*res[1];
    // 	  int myLabel = labels->GetValue(idx) + labelOffsets[processId];

    // 	  if (labels->GetValue(idx) > -1) {
    // 	    if (!uf_find(allLabels, myLabel, neighborLabel)) {
    // 	      uf_unite(allLabels, myLabel, neighborLabel);
    // 	    }
    // 	  }
    // 	}
    // 	++nidx;
    //   }
    // }
    // for (int i = 0; i < numProcesses; ++i) {
    //   recvLengths[i] = numAllLabels;
    //   recvOffsets[i] = i*numAllLabels;
    // }
    // std::vector<int> allLabelUnions(numAllLabels*numProcesses);
    // Controller->AllGatherV(&allLabels[0], &allLabelUnions[0], numAllLabels, 
    // 			   &recvLengths[0], &recvOffsets[0]);

    // for (int i = 0; i < allLabelUnions.size(); ++i) {

    //   int labelId = i%numAllLabels;
    //   if (allLabelUnions[i] != labelId) {
    // 	if (!uf_find(allLabels, allLabelUnions[i], labelId)) {
    // 	  uf_unite(allLabels, allLabelUnions[i], labelId);
    // 	}
    //   }
    // }

    // for (int i = 0; i < allLabels.size(); ++i) {
    //   if (allLabels[i] != i) {
    // 	int rootId = uf_root(allLabels, i);
    // 	allLabels[i] = rootId;
    //   }
    // }
    // std::map<int,int> labelMap;
    // int labelId = 0;
    // for (int i = 0; i < allLabels.size(); ++i) {
    //   if (labelMap.find(allLabels[i]) == labelMap.end()) {
    // 	labelMap[allLabels[i]] = labelId;
    // 	++labelId;
    //   }
    // }

    // float *labels_ptr = labels->GetPointer(0);
    // for (int i = 0; i < labels->GetNumberOfTuples(); ++i) {
    //   if (labels_ptr[i] > -1) {
    // 	int label = labels_ptr[i] + labelOffsets[processId];
    // 	label = allLabels[label];
    // 	labels_ptr[i] = labelMap[label];
    //   }
    // }
  }

  components->SetExtent(vof->GetExtent());
  components->SetXCoordinates(vof->GetXCoordinates());
  components->SetYCoordinates(vof->GetYCoordinates());
  components->SetZCoordinates(vof->GetZCoordinates());
  components->GetCellData()->AddArray(labels);
  components->GetCellData()->SetActiveScalars("Labels");
}

