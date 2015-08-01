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
#include <map>
#include <set>

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
  // set one ghost level -----------------------------------------------------
  int numInputs = this->GetNumberOfInputPorts();
  for (int i = 0; i < numInputs; i++) {
    vtkInformation *inInfo = inputVector[i]->GetInformationObject(0);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);
  }
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);


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

  if (FirstIteration && Controller->GetCommunicator() != 0) {
    // find neighbor processes and global domain bounds
    GetGlobalContext(inInfoVof);
  }

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

      // Stage III -----------------------------------------------------------
      vtkSmartPointer<vtkRectilinearGrid> components = vtkSmartPointer<vtkRectilinearGrid>::New();
      ExtractComponents(inputVof, components);

      // Stage IV ------------------------------------------------------------
      std::vector<float> particleLabels;
      LabelAdvectedParticles(components, particleLabels);

      // Stage V -------------------------------------------------------------
      TransferLabelsToSeeds(particleLabels);

      // Stage VI ------------------------------------------------------------
      vtkSmartPointer<vtkPolyData> boundaries = vtkSmartPointer<vtkPolyData>::New();
      GenerateBoundaries(boundaries);

      // Generate output -----------------------------------------------------
      vtkInformation *outInfo = outputVector->GetInformationObject(0);
      vtkPolyData *output =
      	vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
      output->SetPoints(boundaries->GetPoints());
      output->SetPolys(boundaries->GetPolys());
      output->GetCellData()->AddArray(boundaries->GetCellData()->GetArray("BackLabels"));
      output->GetCellData()->AddArray(boundaries->GetCellData()->GetArray("FrontLabels"));
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
  ParticleIds.clear();
  ParticleProcs.clear();

  Particles.resize(seedPoints->GetNumberOfPoints());
  for (int i = 0; i < seedPoints->GetNumberOfPoints(); ++i) {
    double p[3];
    seedPoints->GetPoint(i, p);
    Particles[i] = make_float4(p[0], p[1], p[2], 1.0f);
  }
  if (Controller->GetCommunicator() != 0) {

    const int processId = Controller->GetLocalProcessId();

    ParticleIds.resize(seedPoints->GetNumberOfPoints());
    ParticleProcs.resize(seedPoints->GetNumberOfPoints());
    for (int i = 0; i < seedPoints->GetNumberOfPoints(); ++i) {
      ParticleIds[i] = i;
      ParticleProcs[i] = processId;
    }
  }
}

//----------------------------------------------------------------------------
void vtkVofTopo::GetGlobalContext(vtkInformation *inInfo)
{
  int numProcesses = Controller->GetNumberOfProcesses();
  std::vector<vtkIdType> RecvLengths(numProcesses);
  std::vector<vtkIdType> RecvOffsets(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    RecvLengths[i] = NUM_SIDES;
    RecvOffsets[i] = i*NUM_SIDES;
  }

  vtkRectilinearGrid *inputVof = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  int LocalExtents[NUM_SIDES];
  inputVof->GetExtent(LocalExtents);

  std::vector<int> AllExtents(NUM_SIDES*numProcesses);
  Controller->AllGatherV(&LocalExtents[0], &AllExtents[0], NUM_SIDES, &RecvLengths[0], &RecvOffsets[0]);

  int GlobalExtents[NUM_SIDES];
  findGlobalExtents(AllExtents, GlobalExtents);

  // reduce extent to one without ghost cells
  if (LocalExtents[0] > GlobalExtents[0])
    LocalExtents[0] += 1;
  if (LocalExtents[1] < GlobalExtents[1])
    LocalExtents[1] -= 1;
  if (LocalExtents[2] > GlobalExtents[2])
    LocalExtents[2] += 1;
  if (LocalExtents[3] < GlobalExtents[3])
    LocalExtents[3] -= 1;
  if (LocalExtents[4] > GlobalExtents[4])
    LocalExtents[4] += 1;
  if (LocalExtents[5] < GlobalExtents[5])
    LocalExtents[5] -= 1;

  // send extents again
  Controller->AllGatherV(&LocalExtents[0], &AllExtents[0], NUM_SIDES, &RecvLengths[0], &RecvOffsets[0]);

  NeighborProcesses.clear();
  NeighborProcesses.resize(NUM_SIDES);
  findNeighbors(LocalExtents, GlobalExtents, AllExtents, NeighborProcesses);

  NumNeighbors = 0;
  for (int i = 0; i < NeighborProcesses.size(); ++i) {
    NumNeighbors += NeighborProcesses[i].size();
  }

  // find domain bounds
  inputVof->GetBounds(&LocalBounds[0]);
  std::vector<double> AllBounds(NUM_SIDES*numProcesses);
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
  std::vector<std::vector<int> > particleIdsToSend(numProcesses);
  std::vector<std::vector<short> > particleProcsToSend(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    particlesToSend[i].resize(0);
    particleIdsToSend[i].resize(0);
    particleProcsToSend[i].resize(0);
  }

  std::vector<float4>::iterator it;
  std::vector<float4> particlesToKeep;
  std::vector<int> particleIdsToKeep;
  std::vector<short> particleProcsToKeep;

  for (int i = 0; i < Particles.size(); ++i) {

    int bound = outOfBounds(Particles[i], LocalBounds, GlobalBounds);
    if (bound > -1) {
      // for (int j = 0; j < NeighborProcesses[bound].size(); ++j) {

      // 	int neighborId = NeighborProcesses[bound][j];
      // 	particlesToSend[neighborId].push_back(Particles[i]);
      // 	particleIdsToSend[neighborId].push_back(ParticleIds[i]);
      // 	particleProcsToSend[neighborId].push_back(ParticleProcs[i]);
      // }
      for (int j = 0; j < numProcesses; ++j) {

  	int neighborId = j;
	if (neighborId != processId) {
	  particlesToSend[neighborId].push_back(Particles[i]);
	  particleIdsToSend[neighborId].push_back(ParticleIds[i]);
	  particleProcsToSend[neighborId].push_back(ParticleProcs[i]);
	}
      }

    }
    else {
      particlesToKeep.push_back(Particles[i]);
      particleIdsToKeep.push_back(ParticleIds[i]);
      particleProcsToKeep.push_back(ParticleProcs[i]);
    }
  }
  Particles = particlesToKeep;
  ParticleIds = particleIdsToKeep;
  ParticleProcs = particleProcsToKeep;

  std::vector<float4> particlesToRecv;
  std::vector<int> particleIdsToRecv;
  std::vector<short> particleProcsToRecv;
  sendData(particlesToSend, particlesToRecv, numProcesses, Controller);
  sendData(particleIdsToSend, particleIdsToRecv, numProcesses, Controller);
  sendData(particleProcsToSend, particleProcsToRecv, numProcesses, Controller);

  // insert the paricles that are within the domain
  for (int i = 0; i < particlesToRecv.size(); ++i) {
    int within = withinBounds(particlesToRecv[i], LocalBounds);
    if (within) {
      Particles.push_back(particlesToRecv[i]);
      ParticleIds.push_back(particleIdsToRecv[i]);
      ParticleProcs.push_back(particleProcsToRecv[i]);
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
      recvLengths[i] = NUM_SIDES;
      recvOffsets[i] = i*NUM_SIDES;
    }

    // -----------------------------------------------------------------------
    // prepare labelled cells to send to neighbors
    int myExtent[NUM_SIDES];
    vof->GetExtent(myExtent);

    std::vector<std::vector<float4> > labelsToSend(6);
    prepareLabelsToSend(NeighborProcesses, myExtent, cellRes, labels, labelsToSend);

    // -----------------------------------------------------------------------
    // send header to neighbors with the number of labels to be send
    int numLabelsToSend[NUM_SIDES];
    for (int i = 0; i < NUM_SIDES; ++i) {
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
    for (int i = 0; i < NUM_SIDES; ++i) {
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
    for (int i = 0; i < NUM_SIDES; ++i) {
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
    for (int i = 0; i < NUM_SIDES; ++i) {
      for (int j = 0; j < NeighborProcesses[i].size(); ++j) {
    	const int RECV_LABEL_DATA_TAG = 100+NeighborProcesses[i][j];

    	Controller->NoBlockReceive((char*)&(labelsToRecv[nidx][0]),
    				   labelsToRecv[nidx].size()*sizeof(float4),
    				   NeighborProcesses[i][j], RECV_LABEL_DATA_TAG, reqs[nidx]);
    	++nidx;
      }
    }
    Controller->WaitAll(NumNeighbors, reqs);

    // -----------------------------------------------------------------------
    // identify equivalent labels from neighbor processes
    std::vector<int> allLabels(numAllLabels);
    for (int i = 0; i < allLabels.size(); ++i) {
      allLabels[i] = i;
    }
    unifyLabelsInProcess(NeighborProcesses, myExtent, cellRes,
			 labels, labelsToRecv, labelOffsets, processId,
			 allLabels);

    for (int i = 0; i < numProcesses; ++i) {
      recvLengths[i] = numAllLabels;
      recvOffsets[i] = i*numAllLabels;
    }
    std::vector<int> allLabelUnions(numAllLabels*numProcesses);
    Controller->AllGatherV(&allLabels[0], &allLabelUnions[0], numAllLabels,
    			   &recvLengths[0], &recvOffsets[0]);

    unifyLabelsInDomain(allLabelUnions, numAllLabels, allLabels, labels,
			labelOffsets, processId);
  }

  components->SetExtent(vof->GetExtent());
  components->SetXCoordinates(vof->GetXCoordinates());
  components->SetYCoordinates(vof->GetYCoordinates());
  components->SetZCoordinates(vof->GetZCoordinates());
  components->GetCellData()->AddArray(labels);
  components->GetCellData()->SetActiveScalars("Labels");
}

//----------------------------------------------------------------------------
void vtkVofTopo::LabelAdvectedParticles(vtkRectilinearGrid *components,
					std::vector<float> &labels)
{
  labels.resize(Particles.size());

  vtkDataArray *data =
    components->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  int nodeRes[3];
  components->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

  for (int i = 0; i < Particles.size(); ++i) {

    if (Particles[i].w == 0.0f) {
      labels[i] = -1.0f;
      continue;
    }
    
    double x[3] = {Particles[i].x, Particles[i].y, Particles[i].z};
    int ijk[3];
    double pcoords[3];
    int particleInsideGrid = components->ComputeStructuredCoordinates(x, ijk, pcoords);

    if (particleInsideGrid) {
      
      int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      float label = data->GetComponent(idx,0);
      labels[i] = label;
    }
    else {
      labels[i] = -1.0f;
    }
  }
}

//----------------------------------------------------------------------------
void vtkVofTopo::TransferLabelsToSeeds(std::vector<float> &particleLabels)
{
  vtkFloatArray *labelsArray = vtkFloatArray::New();
  labelsArray->SetName("Labels");
  labelsArray->SetNumberOfComponents(1);
  labelsArray->SetNumberOfTuples(Seeds->GetNumberOfPoints());
  for (int i = 0; i < Seeds->GetNumberOfPoints(); ++i) {
    labelsArray->SetValue(i, -10.0f);
  }

  if (sizeof(float) != sizeof(int)) {
    vtkDebugMacro("offsets computed assuming same size of int and \
 float, but they have different size");
  }

  if (Controller->GetCommunicator() == 0) {
    for (int i = 0; i < particleLabels.size(); ++i) {
      labelsArray->SetValue(i, particleLabels[i]);
    }
  }
  else {
    
    const int processId = Controller->GetLocalProcessId();
    const int numProcesses = Controller->GetNumberOfProcesses();

    std::vector<std::vector<float> > labelsToSend(numProcesses);
    std::vector<std::vector<int> > idsToSend(numProcesses);
    for (int i = 0; i < numProcesses; ++i) {
      labelsToSend[i].resize(0);
      idsToSend[i].resize(0);
    }

    for (int i = 0; i < particleLabels.size(); ++i) {

      // particle started from a seed in other process - its label and id
      // will be sent to that process
      if (processId != ParticleProcs[i]) {
	labelsToSend[ParticleProcs[i]].push_back(particleLabels[i]);
	idsToSend[ParticleProcs[i]].push_back(ParticleIds[i]);
      }
      // particle started from a seed in this process
      else {
	labelsArray->SetValue(ParticleIds[i], particleLabels[i]);
      }
    }

    // send labels to particle seeds
    std::vector<int> numLabelsToSend(numProcesses);
    std::vector<int> numLabelsToRecv(numProcesses);
    std::vector<float> allLabelsToSend;
    std::vector<int> allIdsToSend;
    allLabelsToSend.resize(0);
    allIdsToSend.resize(0);
    for (int i = 0; i < numProcesses; ++i) {
      numLabelsToSend[i] = labelsToSend[i].size();
      numLabelsToRecv[i] = 0;      
      for (int j = 0; j < labelsToSend[i].size(); ++j) {
	allLabelsToSend.push_back(labelsToSend[i][j]);
	allIdsToSend.push_back(idsToSend[i][j]);
      }
    }

    std::vector<int> RecvLengths(numProcesses);
    std::vector<int> RecvOffsets(numProcesses);
    int numAllLabelsToRecv = 0;
    for (int i = 0; i < numProcesses; ++i) {
      Controller->Scatter((int*)&numLabelsToSend[0], (int*)&numLabelsToRecv[i], 1, i);

      RecvOffsets[i] = numAllLabelsToRecv;
      RecvLengths[i] = numLabelsToRecv[i]*sizeof(float);
      numAllLabelsToRecv += numLabelsToRecv[i];
    }

    std::vector<float> labelsToRecv(numAllLabelsToRecv);
    std::vector<int> idsToRecv(numAllLabelsToRecv, -10000);
    std::vector<vtkIdType> SendLengths(numProcesses);
    std::vector<vtkIdType> SendOffsets(numProcesses);
    int offset = 0;
    for (int i = 0; i < numProcesses; ++i) {
      SendLengths[i] = numLabelsToSend[i]*sizeof(float);
      SendOffsets[i] = offset;
      offset += numLabelsToSend[i]*sizeof(float);
    }

    for (int i = 0; i < numProcesses; ++i) {
      Controller->ScatterV((char*)&allLabelsToSend[0], (char*)&labelsToRecv[RecvOffsets[i]], 
			   &SendLengths[0], &SendOffsets[0], RecvLengths[i], i);
    }
    for (int i = 0; i < numProcesses; ++i) {
      Controller->ScatterV((char*)&allIdsToSend[0], (char*)&idsToRecv[RecvOffsets[i]], 
			   &SendLengths[0], &SendOffsets[0], RecvLengths[i], i);
    }
    for (int i = 0; i < labelsToRecv.size(); ++i) {
      labelsArray->SetValue(idsToRecv[i], labelsToRecv[i]);
    }
  }
  Seeds->GetPointData()->AddArray(labelsArray);
}

//----------------------------------------------------------------------------
void vtkVofTopo::GenerateBoundaries(vtkPolyData *boundaries)
{
  vtkPoints *points = Seeds->GetPoints();
  vtkFloatArray *labels = vtkFloatArray::
    SafeDownCast(Seeds->GetPointData()->GetArray("Labels"));
  vtkIntArray *connectivity = vtkIntArray::
    SafeDownCast(Seeds->GetPointData()->GetArray("Connectivity"));
  vtkShortArray *coords = vtkShortArray::
    SafeDownCast(Seeds->GetPointData()->GetArray("Coords"));

  if (!labels || !connectivity || !coords) {
    vtkDebugMacro("One of the input attributes is not present");
    return;
  }

  generateBoundaries(points, labels, connectivity, coords, boundaries);
}
