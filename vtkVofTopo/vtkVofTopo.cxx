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
#include "vtkMultiBlockDataSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetSurfaceFilter.h"

#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLRectilinearGridWriter.h"1

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <unordered_set>
#include <array>
#include <string>
#include <ctime>
#include <omp.h>

vtkStandardNewMacro(vtkVofTopo);

#include <sys/types.h>
#include <sys/sysinfo.h>

//----------------------------------------------------------------------------
int vtkVofTopo::RequestInformation(vtkInformation *vtkNotUsed(request),
				   vtkInformationVector **inputVector,
				   vtkInformationVector *outputVector)
{
  // optional input port with seeds 
  if (this->GetNumberOfInputConnections(2) > 0) {
    SeedPointsProvided = true;
  }
  
  vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);

  if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {

    unsigned int numberOfInputTimeSteps =
      inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

    this->InputTimeValues.resize(numberOfInputTimeSteps);
    inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
		&this->InputTimeValues[0]);

    if (InputTimeValues.size() > 1 &&
	InputTimeValues[0] > InputTimeValues[1]) {
      Incr = -1.0;
    }
    std::sort(InputTimeValues.begin(), InputTimeValues.end());

    if (numberOfInputTimeSteps == 1) {
      vtkWarningMacro(<<"Not enough input time steps for topology computation");
    }

    if (InitTimeStep < 0) {
      InitTimeStep = 0;
    }
    if (InitTimeStep > InputTimeValues.size()-1) {
      InitTimeStep = InputTimeValues.size()-1;
    }
    if (TargetTimeStep < 0) {
      TargetTimeStep = 0;
    }
    if (TargetTimeStep > InputTimeValues.size()-1) {
      TargetTimeStep = InputTimeValues.size()-1;
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
  const int numInputs = 2;//this->GetNumberOfInputPorts();
  for (int i = 0; i < numInputs; i++) {
    vtkInformation *inInfo = inputVector[i]->GetInformationObject(0);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), NumGhostLevels);
  }
  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), NumGhostLevels);

  if(TimestepT0 == TimestepT1) {

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

    if (LastLoadedTimestep > -1 &&
    	LastLoadedTimestep < TargetTimeStep) {
      UseCache = true;
    }
    else {
      UseCache = false;
      LastLoadedTimestep = -1;
    }

    if (UseCache) {
      ++TimestepT1;
    }
    else {
      TimestepT0 = TimestepT1 = InitTimeStep;
    }
  }
  if (TimestepT1 <= TargetTimeStep) {
    
    int numInputs = 2; //this->GetNumberOfInputPorts();

    for (int i = 0; i < numInputs; i++) {
      vtkInformation *inInfo = inputVector[i]->GetInformationObject(0);

      if (TimestepT1 < static_cast<int>(InputTimeValues.size())) {	
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(),
		    InputTimeValues[TimestepT1]);
	LastLoadedTimestep = TimestepT1;
      }
    }
  }    
  return 1;
}

void writeData(vtkPolyData *data, const int blockId,
	       const int processId, const std::string path)
{
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(data);
  std::string outFile = path;
  outFile += std::to_string(blockId) + std::string("_") + std::to_string(processId);
  outFile += ".vtp";
  writer->SetFileName(outFile.c_str());
  writer->SetDataModeToBinary();
  writer->Write();
}
void writeData(vtkRectilinearGrid *data, const int blockId,
	       const int processId, const std::string path)
{
  vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
  writer->SetInputData(data);
  std::string outFile = path;
  outFile += std::to_string(blockId) + std::string("_") + std::to_string(processId);
  outFile += ".vtr";
  writer->SetFileName(outFile.c_str());
  writer->SetDataModeToBinary();
  writer->Write();
}

//----------------------------------------------------------------------------
int vtkVofTopo::RequestData(vtkInformation *request,
			    vtkInformationVector **inputVector,
			    vtkInformationVector *outputVector)
{
  g_emf0 = EMF0;
  g_emf1 = EMF1;
    
  std::cout << "Timestep T0 T1 = " << TimestepT0 << " " << TimestepT1 << std::endl;
  vtkInformation *inInfoVelocity = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfoVof = inputVector[1]->GetInformationObject(0);
  vtkInformation *inInfoSeeds = nullptr;

  if (SeedPointsProvided) {
    inInfoSeeds = inputVector[2]->GetInformationObject(0);
  }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet *output =
    vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (TimestepT0 == TimestepT1) {
    if (Controller->GetCommunicator() != 0) {
      // find neighbor processes and global domain bounds
      GetGlobalContext(inInfoVof);
    }
    else {
      vtkRectilinearGrid *inputVof = vtkRectilinearGrid::
	SafeDownCast(inInfoVof->Get(vtkDataObject::DATA_OBJECT()));
      // store local extent --------------------------------------------------
      inputVof->GetExtent(LocalExtent);
      inputVof->GetExtent(LocalExtentNoGhosts);
    }
  }

  if (TimestepT0 == TimestepT1) { // first time step
    VofGrid[1]->DeepCopy(vtkRectilinearGrid::
			 SafeDownCast(inInfoVof->Get(vtkDataObject::DATA_OBJECT())));
    VelocityGrid[1]->DeepCopy(vtkRectilinearGrid::
			      SafeDownCast(inInfoVelocity->Get(vtkDataObject::DATA_OBJECT())));
    VofGrid[0]->ShallowCopy(VofGrid[1]);
    VelocityGrid[0]->ShallowCopy(VelocityGrid[1]);
  }
  else {
    VofGrid[0]->ShallowCopy(VofGrid[1]);
    VelocityGrid[0]->ShallowCopy(VelocityGrid[1]);
    
    VofGrid[1]->DeepCopy(vtkRectilinearGrid::
			 SafeDownCast(inInfoVof->Get(vtkDataObject::DATA_OBJECT())));
    VelocityGrid[1]->DeepCopy(vtkRectilinearGrid::
			      SafeDownCast(inInfoVelocity->Get(vtkDataObject::DATA_OBJECT())));
  }
  // Stage I ---------------------------------------------------------------
  if (TimestepT0 == TimestepT1) {
    if (!UseCache) {

      if (SeedPointsProvided) {
	vtkPolyData *seeds = vtkPolyData::
	  SafeDownCast(inInfoSeeds->Get(vtkDataObject::DATA_OBJECT()));
	InitParticles(VofGrid[0], seeds);	
      }
      else {
	InitParticles(VofGrid[0], nullptr);	
      }

      Uncertainty.clear();
      Uncertainty.resize(Particles.size(),0.0f);

      Boundaries->SetPoints(vtkPoints::New());
      vtkCellArray *cells = vtkCellArray::New();
      Boundaries->SetPolys(cells);
    }
  }

  // Stage II --------------------------------------------------------------  
  if (TimestepT0 != TimestepT1) {    
    if(TimestepT0 < TargetTimeStep) {

      // AdvectParticlesInt(VofGrid, VelocityGrid);
      AdvectParticles(VofGrid, VelocityGrid);
      
      if (StoreIntermParticles && TimestepT0 < TargetTimeStep-1) {
	IntermParticles.push_back(Particles);
	IntermParticlesTimeStamps.push_back(TimestepT1);
	if (Controller->GetCommunicator() != 0) {
	  IntermParticleIds.push_back(ParticleIds);
	  IntermParticleProcs.push_back(ParticleProcs);
	}
      }
    }
    
    if (ComputeComponentLabels) {
      bool finishedAdvection = TimestepT1 >= TargetTimeStep;
      if (finishedAdvection) {
	// Stage III -----------------------------------------------------------
	vtkSmartPointer<vtkRectilinearGrid> components = vtkSmartPointer<vtkRectilinearGrid>::New();
	ExtractComponents(VofGrid[1], components);
	
	// Stage IV ------------------------------------------------------------
	std::vector<float> particleLabels;
	LabelAdvectedParticles(components, particleLabels);

	// Stage V -------------------------------------------------------------
	TransferParticleDataToSeeds(particleLabels, "Labels");
	TransferParticleDataToSeeds(Uncertainty, "Uncertainty");

	if (StoreIntermParticles) {
	  TransferIntermParticlesToSeeds(IntermParticles,
					 IntermParticleIds,
					 IntermParticleProcs);
	}
	
	// Transfer seed points from neighbors ---------------------------------
	vtkPolyData *boundarySeeds = vtkPolyData::New();
	if (Controller->GetCommunicator() != 0) {
	  ExchangeBoundarySeedPoints(boundarySeeds);
	}
	
	// Stage VI ------------------------------------------------------------
	GenerateBoundaries(Boundaries, boundarySeeds, this->Seeds);

	boundarySeeds->Delete();

	// Generate output -----------------------------------------------------
	vtkPolyData *particles = vtkPolyData::New();
	vtkPoints *ppoints = vtkPoints::New();
	ppoints->SetNumberOfPoints(Particles.size());
	vtkFloatArray *labels = vtkFloatArray::New();
	labels->SetName("Labels");
	labels->SetNumberOfComponents(1);
	labels->SetNumberOfTuples(particleLabels.size());
	vtkFloatArray *uncertainty = vtkFloatArray::New();
	uncertainty->SetName("Uncertainty");
	uncertainty->SetNumberOfComponents(1);
	uncertainty->SetNumberOfTuples(Uncertainty.size());

	for (int i = 0; i < Particles.size(); ++i) {

	  float p[3] = {Particles[i].x, Particles[i].y, Particles[i].z};
	  ppoints->SetPoint(i, p);
	  labels->SetValue(i, particleLabels[i]);
	  uncertainty->SetValue(i, Uncertainty[i]);
	}
	particles->SetPoints(ppoints);
	particles->GetPointData()->AddArray(labels);
	particles->GetPointData()->AddArray(uncertainty);

	// //
	// vtkFloatArray *ids = vtkFloatArray::New();
	// ids->SetName("Ids");
	// ids->SetNumberOfComponents(1);
	// ids->SetNumberOfTuples(Particles.size());
	// for (int i = 0; i < Particles.size(); ++i) {

	//   ids->SetValue(i, i);
	// }

	// Seeds->GetPointData()->AddArray(ids);
	// particles->GetPointData()->AddArray(ids);
	// //
	
	output->SetBlock(0, Seeds);
	output->SetBlock(1, particles);
	output->SetBlock(2, Boundaries);
	output->SetBlock(3, components);

	// writeData(Seeds, 0, Controller->GetLocalProcessId(), "out2/out2_");
	// writeData(particles, 1, Controller->GetLocalProcessId(), "out2/out2_");
	// writeData(Boundaries, 2, Controller->GetLocalProcessId(), "out2/out2_");
	// writeData(components, 3, Controller->GetLocalProcessId(), "out2/out2_");

	if (StoreIntermParticles) { 
	  int numPoints = 0;
	  for (const auto &ps : IntermParticles) {
	    numPoints += ps.size();
	  }

	  vtkPolyData *intParticles = vtkPolyData::New();
	  vtkPoints *intPoints = vtkPoints::New();
	  intPoints->SetNumberOfPoints(numPoints);
	  
	  vtkFloatArray *intTimeStamps = vtkFloatArray::New();
	  intTimeStamps->SetName("IntermediateTimeStamps");
	  intTimeStamps->SetNumberOfComponents(1);
	  intTimeStamps->SetNumberOfTuples(numPoints);
	  
	  vtkFloatArray *intLabels = vtkFloatArray::New();
	  intLabels->SetName("Labels");
	  intLabels->SetNumberOfComponents(1);
	  intLabels->SetNumberOfTuples(numPoints);
	  
	  int idx = 0;


	  if (Controller->GetCommunicator() == 0) { 
	    for (int i = 0; i < IntermParticles.size(); ++i) {
	      for (int j = 0; j < IntermParticles[i].size(); ++j) {
	      
		const float4 &p = IntermParticles[i][j];
		float pf[3] = {p.x, p.y, p.z};
		intPoints->SetPoint(idx, pf);
		intTimeStamps->SetValue(idx, IntermParticlesTimeStamps[i]);
	      
		float label = Seeds->GetPointData()->GetArray("Labels")->GetComponent(j,0);
		intLabels->SetValue(idx, label);

		++idx;
	      }
	    }
	  }
	  else {
	    for (int i = 0; i < IntermParticles.size(); ++i) {
	      for (int j = 0; j < IntermParticles[i].size(); ++j) {
	      
		const float4 &p = IntermParticles[i][j];
		float pf[3] = {p.x, p.y, p.z};
		intPoints->SetPoint(idx, pf);
		intTimeStamps->SetValue(idx, IntermParticlesTimeStamps[i]);
	      
		const int id = IntermParticleIds[i][j];
		float label = Seeds->GetPointData()->GetArray("Labels")->GetComponent(id,0);
		intLabels->SetValue(idx, label);

		++idx;
	      }
	    }
	  }
	  intParticles->SetPoints(intPoints);
	  intParticles->GetPointData()->AddArray(intTimeStamps);
	  intParticles->GetPointData()->AddArray(intLabels);

	  output->SetBlock(4, intParticles);

	}
	// {
	//   std::vector<float4> vertices;
	//   std::vector<float4> normals;
	//   std::vector<int> indices;
	//   std::vector<float> labels;
	//   std::vector<float> timestamps;
	//   std::vector<std::vector<int>> labelOffsets(IntermParticles.size());

	//   int vertexID = 0;
	//   for (int i = 0; i < IntermParticles.size(); ++i) {
	  
	//     const std::vector<float4> &ip = IntermParticles[i];

	//     generateBoundary(ip, particleLabels, VofGrid[0], Refinement,
	// LocalExtentNoGhosts,
	// 		     vertexID, labelOffsets[i], vertices, normals, indices);
	//   }
	//   {

	//     vtkPoints *outputPoints = vtkPoints::New();
	//     outputPoints->SetNumberOfPoints(vertices.size());
	//     for (int i = 0; i < vertices.size(); ++i) {
    
	//       double p[3] = {vertices[i].x,
	// 		     vertices[i].y,
	// 		     vertices[i].z};
	//       outputPoints->SetPoint(i, p);        
	//     }

	//     vtkIdTypeArray *cells = vtkIdTypeArray::New();
	//     cells->SetNumberOfComponents(1);
	//     cells->SetNumberOfTuples(indices.size()/3*4);
	//     for (int i = 0; i < indices.size()/3; ++i) {
	//       cells->SetValue(i*4+0,3);
	//       cells->SetValue(i*4+1,indices[i*3+0]);
	//       cells->SetValue(i*4+2,indices[i*3+1]);
	//       cells->SetValue(i*4+3,indices[i*3+2]);
	//     }

	//     vtkCellArray *outputTriangles = vtkCellArray::New();
	//     outputTriangles->SetNumberOfCells(indices.size()/3);
	//     outputTriangles->SetCells(indices.size()/3, cells);

	//     // vtkShortArray *boundaryLabels = vtkShortArray::New();
	//     // boundaryLabels->SetName("Labels");
	//     // boundaryLabels->SetNumberOfComponents(1);
	//     // boundaryLabels->SetNumberOfTuples(outputPoints->GetNumberOfPoints());

	//     // std::cout << "labelOffsets.size() = " << labelOffsets.size() << std::endl;
	//     // for (int i = 0; i < labelOffsets.size(); ++i) { // num time steps
	//     //   std::cout << "labelOffsets[" << i << "].size() = " << labelOffsets[i].size() << std::endl;
	//     //   for (int j = 0; j < labelOffsets[i].size()-1; ++j) { // num unique labels
	//     // 	std::cout << labelOffsets[i][j] << " - " << labelOffsets[i][j+1] << std::endl;
	//     // 	for (int k = labelOffsets[i][j]; k < labelOffsets[i][j+1]; ++k) {
	//     // 	  boundaryLabels->SetValue(k, j);
	//     // 	}
	//     //   }
	//     // }

	//     vtkFloatArray *pointNormals = vtkFloatArray::New();
	//     pointNormals->SetName("Normals");
	//     pointNormals->SetNumberOfComponents(3);
	//     pointNormals->SetNumberOfTuples(normals.size());
	//     for (int i = 0; i < normals.size(); ++i) {
	//       float n[3] = {normals[i].x, normals[i].y, normals[i].z};
	//       pointNormals->SetTuple3(i,n[0],n[1],n[2]);
	//     }
	//     vtkPolyData *boundaries = vtkPolyData::New();
	//     boundaries->SetPoints(outputPoints);
	//     boundaries->SetPolys(outputTriangles);
	//     // boundaries->GetPointData()->AddArray(boundaryLabels);
	//     boundaries->GetPointData()->SetNormals(pointNormals);

	//     output->SetBlock(4, boundaries);
	//   }
	// }
      }
    }
  }

  TimestepT0 = TimestepT1;      
  bool finishedAdvection = TimestepT1 >= TargetTimeStep;
  if (finishedAdvection) {
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
  }
  else {
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
    ++TimestepT1;
  }
  return 1;
}

//----------------------------------------------------------------------------
void vtkVofTopo::InitParticles(vtkRectilinearGrid *vof, vtkPolyData *seeds)
{
  vtkPoints *seedPoints;

  if (seeds) {
    if (Controller->GetCommunicator() == 0) {
      seedPoints = seeds->GetPoints();
    }
    else {

      seedPoints = vtkPoints::New();
      
      // send seed points to all other processes
      Controller->Broadcast(seeds, 0);
      const int numPoints = seeds->GetNumberOfPoints();

      // if running in parallel, get only points in the subdomain
      for (int i = 0; i < numPoints; ++i) {
	double p[3];
	seeds->GetPoint(i, p);
	if (withinBounds(make_float4(p[0],p[1],p[2],0.0f), BoundsNoGhosts)) {
	  seedPoints->InsertNextPoint(p);
	}
      }
    }
  }
  else {
    seedPoints = vtkPoints::New();
    generateSeedPointsPLIC(vof, Refinement, seedPoints, GlobalExtent, NumGhostLevels);
  }
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

  if (Seeds != 0) {
    Seeds->Delete();
  }
  Seeds = vtkPolyData::New();
  Seeds->SetPoints(seedPoints);
}

//----------------------------------------------------------------------------
void vtkVofTopo::GetGlobalContext(vtkInformation *inInfo)
{
  int processId = Controller->GetLocalProcessId();
  int numProcesses = Controller->GetNumberOfProcesses();

  // prepare buffers for communication ---------------------------------------
  std::vector<vtkIdType> RecvLengths(numProcesses);
  std::vector<vtkIdType> RecvOffsets(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    RecvLengths[i] = NUM_SIDES;
    RecvOffsets[i] = i*NUM_SIDES;
  }

  vtkRectilinearGrid *inputVof = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // store local extent ------------------------------------------------------
  inputVof->GetExtent(LocalExtent);
  
  // find global extent ------------------------------------------------------
  std::vector<int> AllExtents(NUM_SIDES*numProcesses);
  Controller->AllGatherV(&LocalExtent[0], &AllExtents[0], NUM_SIDES, 
			 &RecvLengths[0], &RecvOffsets[0]);
  findGlobalExtent(AllExtents, GlobalExtent);

  // reduce extent to one without ghost cells --------------------------------
  inputVof->GetExtent(LocalExtentNoGhosts);
  if (LocalExtent[0] > GlobalExtent[0]) LocalExtentNoGhosts[0] += NumGhostLevels;
  if (LocalExtent[1] < GlobalExtent[1]) LocalExtentNoGhosts[1] -= NumGhostLevels;
  if (LocalExtent[2] > GlobalExtent[2]) LocalExtentNoGhosts[2] += NumGhostLevels;
  if (LocalExtent[3] < GlobalExtent[3]) LocalExtentNoGhosts[3] -= NumGhostLevels;
  if (LocalExtent[4] > GlobalExtent[4]) LocalExtentNoGhosts[4] += NumGhostLevels;
  if (LocalExtent[5] < GlobalExtent[5]) LocalExtentNoGhosts[5] -= NumGhostLevels;

  // find neighboring subdomains ---------------------------------------------
  std::vector<int> AllExtentsNoGhosts(NUM_SIDES*numProcesses);
  Controller->AllGatherV(&LocalExtentNoGhosts[0], &AllExtentsNoGhosts[0], 
			 NUM_SIDES, &RecvLengths[0], &RecvOffsets[0]);
  NeighborProcesses.clear();
  NeighborProcesses.resize(NUM_SIDES);
  findNeighbors(LocalExtentNoGhosts, GlobalExtent, AllExtentsNoGhosts, 
		NeighborProcesses, processId);  

  NumNeighbors = 0;
  for (int i = 0; i < NeighborProcesses.size(); ++i) {
    NumNeighbors += NeighborProcesses[i].size();
  }

  // find domain bounds ------------------------------------------------------
  // local bounds ------------------------------------------------------------
  inputVof->GetBounds(&LocalBounds[0]);

  // global bounds -----------------------------------------------------------
  std::vector<double> AllBounds(NUM_SIDES*numProcesses);
  Controller->AllGatherV(&LocalBounds[0], &AllBounds[0], 6, &RecvLengths[0], &RecvOffsets[0]);
  findGlobalBounds(AllBounds, GlobalBounds);

  // local bounds without ghosts ---------------------------------------------
  vtkDataArray *coords[3] = {inputVof->GetXCoordinates(), 
			     inputVof->GetYCoordinates(), 
			     inputVof->GetZCoordinates()};
  int ncrds[3] = {coords[0]->GetNumberOfTuples(),
		      coords[1]->GetNumberOfTuples(),
		      coords[2]->GetNumberOfTuples()};

  BoundsNoGhosts[0] = coords[0]->GetComponent(LocalExtentNoGhosts[0]-LocalExtent[0], 0);
  BoundsNoGhosts[1]=coords[0]->GetComponent(ncrds[0]-1-(LocalExtent[1]-LocalExtentNoGhosts[1]), 0);
  BoundsNoGhosts[2] = coords[1]->GetComponent(LocalExtentNoGhosts[2]-LocalExtent[2], 0);
  BoundsNoGhosts[3]=coords[1]->GetComponent(ncrds[1]-1-(LocalExtent[3]-LocalExtentNoGhosts[3]), 0);
  BoundsNoGhosts[4] = coords[2]->GetComponent(LocalExtentNoGhosts[4]-LocalExtent[4], 0);
  BoundsNoGhosts[5]=coords[2]->GetComponent(ncrds[2]-1-(LocalExtent[5]-LocalExtentNoGhosts[5]), 0);
}

//----------------------------------------------------------------------------
void vtkVofTopo::AdvectParticles(vtkRectilinearGrid *vof[2],
				 vtkRectilinearGrid *velocity[2])
{
  float dt = InputTimeValues[TimestepT1] - InputTimeValues[TimestepT0];
  if (TimeStepDelta != 0.0) {
    dt = TimeStepDelta;
  }
  
  dt *= Incr;

  advectParticles(vof, velocity, Particles, Uncertainty, dt,
		  IntegrationMethod, PLICCorrection, VOFCorrection, RK4NumSteps);
  if (Controller->GetCommunicator() != 0) {
    ExchangeParticles();
  }
}

//----------------------------------------------------------------------------
void vtkVofTopo::AdvectParticlesInt(vtkRectilinearGrid *vof[2],
				    vtkRectilinearGrid *velocity[2])
{
  float dt = InputTimeValues[TimestepT1] - InputTimeValues[TimestepT0];
  if (TimeStepDelta != 0.0) {
    dt = TimeStepDelta;
  }
  
  dt *= Incr;

  vtkSmartPointer<vtkRectilinearGrid> intVof = vtkSmartPointer<vtkRectilinearGrid>::New();
  vtkSmartPointer<vtkRectilinearGrid> intVelocity = vtkSmartPointer<vtkRectilinearGrid>::New();
      
  InterpolateField(VofGrid, VelocityGrid, intVof, intVelocity, 0.5f);

  {
  vtkRectilinearGrid *vofInt[3] = {vof[0], intVof, vof[1]};
  vtkRectilinearGrid *velocityInt[3] = {velocity[0], intVelocity, velocity[1]};

  advectParticlesInt(vofInt, velocityInt, Particles, Uncertainty, dt, PLICCorrection, VOFCorrection);
  if (Controller->GetCommunicator() != 0) {
    ExchangeParticles();
  }
  }
  
  // {
  //   vtkRectilinearGrid *vofTmp[2] = {vof[0], intVof};
  //   vtkRectilinearGrid *velocityTmp[2] = {velocity[0], intVelocity};
  //   advectParticles(vofTmp, velocityTmp, Particles, Uncertainty, dt/2.0f,
  // 		    IntegrationMethod, PLICCorrection, VOFCorrection, RK4NumSteps);
  //   if (Controller->GetCommunicator() != 0) {
  //     ExchangeParticles();
  //   }
  // }  
  // {
  //   vtkRectilinearGrid *vofTmp[2] = {intVof,vof[1]};
  //   vtkRectilinearGrid *velocityTmp[2] = {intVelocity,velocity[1]};
  //   advectParticles(vofTmp, velocityTmp, Particles, Uncertainty, dt/2.0f,
  // 		    IntegrationMethod, PLICCorrection, VOFCorrection, RK4NumSteps);
  //   if (Controller->GetCommunicator() != 0) {
  //     ExchangeParticles();
  //   }
  // }  
}

//----------------------------------------------------------------------------
void vtkVofTopo::ExchangeParticles()
{
  int numProcesses = Controller->GetNumberOfProcesses();
  int processId = Controller->GetLocalProcessId();

  // one vector for each side of the process
  std::vector<std::vector<float4>> particlesToSend(numProcesses);
  std::vector<std::vector<int> > particleIdsToSend(numProcesses);
  std::vector<std::vector<short> > particleProcsToSend(numProcesses);
  std::vector<std::vector<float>> uncertaintyToSend(numProcesses);
  
  for (int i = 0; i < numProcesses; ++i) {
    particlesToSend[i].resize(0);
    particleIdsToSend[i].resize(0);
    particleProcsToSend[i].resize(0);
    uncertaintyToSend[i].resize(0);
  }

  std::vector<float4> particlesToKeep;
  std::vector<int> particleIdsToKeep;
  std::vector<short> particleProcsToKeep;
  std::vector<float> uncertaintyToKeep;

  for (int i = 0; i < Particles.size(); ++i) {

    int bound = outOfBounds(Particles[i], BoundsNoGhosts, GlobalBounds);
    if (bound > -1) {
      for (int j = 0; j < NeighborProcesses[bound].size(); ++j) {
      // for (int j = 0; j < numProcesses; ++j) {	

  	int neighborId = NeighborProcesses[bound][j];
	if (neighborId != processId) {
	  particlesToSend[neighborId].push_back(Particles[i]);
	  particleIdsToSend[neighborId].push_back(ParticleIds[i]);
	  particleProcsToSend[neighborId].push_back(ParticleProcs[i]);
	  uncertaintyToSend[neighborId].push_back(Uncertainty[i]);
	}
      }
    }
    else {
      particlesToKeep.push_back(Particles[i]);
      particleIdsToKeep.push_back(ParticleIds[i]);
      particleProcsToKeep.push_back(ParticleProcs[i]);
      uncertaintyToKeep.push_back(Uncertainty[i]);
    }
  }
 
  Particles = particlesToKeep;
  ParticleIds = particleIdsToKeep;
  ParticleProcs = particleProcsToKeep;
  Uncertainty = uncertaintyToKeep;

  std::vector<float4> particlesToRecv;
  std::vector<int> particleIdsToRecv;
  std::vector<short> particleProcsToRecv;
  std::vector<float> uncertaintyToRecv;
  sendData(particlesToSend, particlesToRecv, numProcesses, Controller);
  sendData(particleIdsToSend, particleIdsToRecv, numProcesses, Controller);
  sendData(particleProcsToSend, particleProcsToRecv, numProcesses, Controller);
  sendData(uncertaintyToSend, uncertaintyToRecv, numProcesses, Controller);

  // insert the paricles that are within the domain
  for (int i = 0; i < particlesToRecv.size(); ++i) {
    int within = withinBounds(particlesToRecv[i], BoundsNoGhosts);
    if (within) {
      Particles.push_back(particlesToRecv[i]);
      ParticleIds.push_back(particleIdsToRecv[i]);
      ParticleProcs.push_back(particleProcsToRecv[i]);
      Uncertainty.push_back(uncertaintyToRecv[i]);
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
    prepareLabelsToSend(NeighborProcesses, myExtent, GlobalExtent,
			cellRes, labels, labelsToSend, NumGhostLevels);

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

  components->CopyStructure(vof);
  components->GetCellData()->AddArray(labels);
  components->GetCellData()->SetActiveScalars("Labels");

  int myExtent[NUM_SIDES];
  vof->GetExtent(myExtent);
  int updatedExtent[6] = {myExtent[0]+NumGhostLevels,
			  myExtent[1]-NumGhostLevels,
			  myExtent[2]+NumGhostLevels,
			  myExtent[3]-NumGhostLevels,
			  myExtent[4]+NumGhostLevels,
			  myExtent[5]-NumGhostLevels};
  components->Crop(updatedExtent);
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

    if (Particles[i].w <= g_emf0) {
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
void vtkVofTopo::TransferParticleDataToSeeds(std::vector<float> &particleData,
					     const std::string arrayName)
{
  vtkFloatArray *dataArray = vtkFloatArray::New();
  dataArray->SetName(arrayName.c_str());
  dataArray->SetNumberOfComponents(1);
  dataArray->SetNumberOfTuples(Seeds->GetNumberOfPoints());
  for (int i = 0; i < Seeds->GetNumberOfPoints(); ++i) {
    dataArray->SetValue(i, -10.0f);
  }

  if (sizeof(float) != sizeof(int)) {
    vtkDebugMacro("offsets computed assuming same size of int and \
 float, but they have different size");
  }

  if (Controller->GetCommunicator() == 0) {
    for (int i = 0; i < particleData.size(); ++i) {
      dataArray->SetValue(i, particleData[i]);
    }
  }
  else { // parallel
    
    const int processId = Controller->GetLocalProcessId();
    const int numProcesses = Controller->GetNumberOfProcesses();

    std::vector<std::vector<float>> dataToSend(numProcesses);
    std::vector<std::vector<int>> idsToSend(numProcesses);
    for (int i = 0; i < numProcesses; ++i) {
      dataToSend[i].resize(0);
      idsToSend[i].resize(0);
    }

    for (int i = 0; i < particleData.size(); ++i) {

      // particle started from a seed in other process - its label and id
      // will be sent to that process
      if (processId != ParticleProcs[i]) {
	dataToSend[ParticleProcs[i]].push_back(particleData[i]);
	idsToSend[ParticleProcs[i]].push_back(ParticleIds[i]);
      }
      // particle started from a seed in this process
      else {
	dataArray->SetValue(ParticleIds[i], particleData[i]);
      }
    }

    // send data to particle seeds
    std::vector<int> numDataToSend(numProcesses);
    std::vector<int> numDataToRecv(numProcesses);
    std::vector<float> allDataToSend;
    std::vector<int> allIdsToSend;
    allDataToSend.resize(0);
    allIdsToSend.resize(0);
    for (int i = 0; i < numProcesses; ++i) {
      numDataToSend[i] = dataToSend[i].size();
      numDataToRecv[i] = 0;      
      for (int j = 0; j < dataToSend[i].size(); ++j) {
	allDataToSend.push_back(dataToSend[i][j]);
	allIdsToSend.push_back(idsToSend[i][j]);
      }
    }

    // no longer needed
    dataToSend.clear();
    idsToSend.clear();

    std::vector<int> RecvLengths(numProcesses);
    std::vector<int> RecvOffsets(numProcesses);
    int numAllDataToRecv = 0;
    for (int i = 0; i < numProcesses; ++i) {
      Controller->Scatter((int*)&numDataToSend[0], (int*)&numDataToRecv[i], 1, i);

      RecvOffsets[i] = numAllDataToRecv;
      RecvLengths[i] = numDataToRecv[i]*sizeof(float);
      numAllDataToRecv += numDataToRecv[i];
    }

    std::vector<float> dataToRecv(numAllDataToRecv);
    std::vector<int> idsToRecv(numAllDataToRecv, -10000);
    std::vector<vtkIdType> SendLengths(numProcesses);
    std::vector<vtkIdType> SendOffsets(numProcesses);
    int offset = 0;
    for (int i = 0; i < numProcesses; ++i) {
      SendLengths[i] = numDataToSend[i]*sizeof(float);
      SendOffsets[i] = offset;
      offset += numDataToSend[i]*sizeof(float);
    }

    for (int i = 0; i < numProcesses; ++i) {
      Controller->ScatterV((char*)&allDataToSend[0], (char*)&dataToRecv[RecvOffsets[i]], 
			   &SendLengths[0], &SendOffsets[0], RecvLengths[i], i);
    }
    for (int i = 0; i < numProcesses; ++i) {
      Controller->ScatterV((char*)&allIdsToSend[0], (char*)&idsToRecv[RecvOffsets[i]], 
			   &SendLengths[0], &SendOffsets[0], RecvLengths[i], i);
    }
    for (int i = 0; i < dataToRecv.size(); ++i) {
      dataArray->SetValue(idsToRecv[i], dataToRecv[i]);
    }

  }
  Seeds->GetPointData()->AddArray(dataArray);
}

//----------------------------------------------------------------------------
void vtkVofTopo::GenerateBoundaries(vtkPolyData *boundaries, vtkPolyData *boundarySeeds,
				    vtkPolyData *seeds)
{
  vtkPoints *points = seeds->GetPoints();
  vtkFloatArray *labels = vtkFloatArray::
    SafeDownCast(seeds->GetPointData()->GetArray("Labels"));

  vtkPoints *boundarySeedPoints = 0;
  vtkFloatArray *boundarySeedLabels = 0;
  if (boundarySeeds != 0) {
    boundarySeedPoints = boundarySeeds->GetPoints();
    boundarySeedLabels = vtkFloatArray::SafeDownCast(boundarySeeds->GetPointData()->GetArray("Labels"));
  }
  generateBoundaries(points, boundarySeedPoints,
		     labels, boundarySeedLabels,
		     this->VofGrid[1], boundaries, 
		     this->LocalExtentNoGhosts,
		     this->Refinement);
}

int isInside(const int ijk[3], const int cellRes[3],
	     const int numGhostLevels, const int outerBoundary)
{
  if (ijk[0] >= numGhostLevels-outerBoundary &&
      ijk[0] < cellRes[0]-numGhostLevels+outerBoundary &&
      ijk[1] >= numGhostLevels-outerBoundary &&
      ijk[1] < cellRes[1]-numGhostLevels+outerBoundary &&
      ijk[2] >= numGhostLevels-outerBoundary &&
      ijk[2] < cellRes[2]-numGhostLevels+outerBoundary) {
    return 1;
  } 
  return 0;
}

//----------------------------------------------------------------------------
void vtkVofTopo::ExchangeBoundarySeedPoints(vtkPolyData *boundarySeeds)
{
  const int boundarySize = 2;
  int extent[6];
  this->VofGrid[1]->GetExtent(extent);

  int nodeRes[3] = {this->VofGrid[1]->GetXCoordinates()->GetNumberOfTuples(), 
		    this->VofGrid[1]->GetYCoordinates()->GetNumberOfTuples(), 
		    this->VofGrid[1]->GetZCoordinates()->GetNumberOfTuples()};
  int cellRes[3] = {nodeRes[0]-1,nodeRes[1]-1,nodeRes[2]-1};
  
  // indices of inner cells (without domain offset)
  int innerExtent[6] = {extent[0] > GlobalExtent[0] ? NumGhostLevels : 0,
			extent[1] < GlobalExtent[1] ? cellRes[0]-1-NumGhostLevels : cellRes[0]-1,
			extent[2] > GlobalExtent[2] ? NumGhostLevels : 0,
			extent[3] < GlobalExtent[3] ? cellRes[1]-1-NumGhostLevels : cellRes[1]-1,
			extent[4] > GlobalExtent[4] ? NumGhostLevels : 0,
			extent[5] < GlobalExtent[5] ? cellRes[2]-1-NumGhostLevels : cellRes[2]-1};

  int bids[6] = {extent[0] > GlobalExtent[0] ? innerExtent[0] : -1,
		 extent[1] < GlobalExtent[1] ? innerExtent[1] : -1,
		 extent[2] > GlobalExtent[2] ? innerExtent[2] : -1,
		 extent[3] < GlobalExtent[3] ? innerExtent[3] : -1,
		 extent[4] > GlobalExtent[4] ? innerExtent[4] : -1,
		 extent[5] < GlobalExtent[5] ? innerExtent[5] : -1};

  std::set<std::array<int,3>> boundaryCells;
  boundaryCells.clear();

  if (bids[0] > -1) { // left
    for (int i = bids[0]; i < bids[0]+boundarySize; ++i)
      for (int k = innerExtent[4]; k <= innerExtent[5]; ++k)
	for (int j = innerExtent[2]; j <= innerExtent[3]; ++j)
	  boundaryCells.emplace(std::array<int,3>{i,j,k});
  }
  if (bids[1] > -1) { // right
    for (int i = bids[1]; i > bids[1]-boundarySize; --i) 
      for (int k = innerExtent[4]; k <= innerExtent[5]; ++k)
	for (int j = innerExtent[2]; j <= innerExtent[3]; ++j)
	  boundaryCells.emplace(std::array<int,3>{i,j,k});
  }
  if (bids[2] > -1) { // bottom
    for (int j = bids[2]; j < bids[2]+boundarySize; ++j)
      for (int k = innerExtent[4]; k <= innerExtent[5]; ++k)
	for (int i = innerExtent[0]; i <= innerExtent[1]; ++i)
	  boundaryCells.emplace(std::array<int,3>{i,j,k});
  }
  if (bids[3] > -1) { // top
    for (int j = bids[3]; j > bids[3]-boundarySize; --j)
      for (int k = innerExtent[4]; k <= innerExtent[5]; ++k)
	for (int i = innerExtent[0]; i <= innerExtent[1]; ++i)
	  boundaryCells.emplace(std::array<int,3>{i,j,k});
  }
  if (bids[4] > -1) { // back
    for (int k = bids[4]; k < bids[4]+boundarySize; ++k)
      for (int j = innerExtent[2]; j <= innerExtent[3]; ++j)
	for (int i = innerExtent[0]; i <= innerExtent[1]; ++i)
	  boundaryCells.emplace(std::array<int,3>{i,j,k});
  }
  if (bids[5] > -1) { // front
    for (int k = bids[5]; k > bids[5]-boundarySize; --k)
      for (int j = innerExtent[2]; j <= innerExtent[3]; ++j)
	for (int i = innerExtent[0]; i <= innerExtent[1]; ++i)
	  boundaryCells.emplace(std::array<int,3>{i,j,k});
  }

  vtkPoints *seedPoints = Seeds->GetPoints();
  vtkFloatArray *labels = vtkFloatArray::SafeDownCast(Seeds->GetPointData()->GetArray("Labels"));
  const int numSeedPoints = seedPoints->GetNumberOfPoints();

  std::vector<float3> pointsToSend;
  std::vector<float> labelsToSend;
  pointsToSend.clear();
  labelsToSend.clear();
  
  for (int i = 0; i < numSeedPoints; ++i) {
    int ijk[3];
    double pcoords[3];
    double x[3];
    seedPoints->GetPoint(i, x);
    VofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);

    if (boundaryCells.find(std::array<int,3>{ijk[0],ijk[1],ijk[2]}) != boundaryCells.end()) {
      pointsToSend.push_back(make_float3(x[0],x[1],x[2]));
      labelsToSend.push_back(labels->GetValue(i));
    }
  }

  int numProcesses = Controller->GetNumberOfProcesses();
  std::vector<float3> pointsToRecv;
  std::vector<float> labelsToRecv;
  pointsToRecv.clear();
  labelsToRecv.clear();
  sendData(pointsToSend, pointsToRecv, numProcesses, this->Controller);
  sendData(labelsToSend, labelsToRecv, numProcesses, this->Controller);

  vtkPoints *boundarySeedPoints = vtkPoints::New();
  vtkFloatArray *boundarySeedLabels = vtkFloatArray::New();
  boundarySeedLabels->SetName("Labels");
  boundarySeedLabels->SetNumberOfComponents(1);
  
  for (int i = 0; i < pointsToRecv.size(); ++i) {
    double x[3] = {pointsToRecv[i].x, pointsToRecv[i].y, pointsToRecv[i].z};
    int ijk[3];
    double pcoords[3];
    int inside = VofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
    int insideGrid = isInside(ijk,cellRes,NumGhostLevels,boundarySize);
    if (inside && insideGrid) {
    
      boundarySeedPoints->InsertNextPoint(x);
      boundarySeedLabels->InsertNextTuple1(labelsToRecv[i]);
    }
  }

  boundarySeeds->SetPoints(boundarySeedPoints);
  boundarySeeds->GetPointData()->AddArray(boundarySeedLabels);
}



//----------------------------------------------------------------------------
vtkVofTopo::vtkVofTopo() :
  LastLoadedTimestep(-1),
  UseCache(false),
  IterType(ITERATE_OVER_TARGET),
  ComputeComponentLabels(1),
  Seeds(0),
  Incr(1.0),
  TimestepT0(-1),
  TimestepT1(-1),
  NumGhostLevels(4),
  SeedPointsProvided(false),
  IntegrationMethod(0), // Heun
  PLICCorrection(0),
  VOFCorrection(0),
  RK4NumSteps(8),
  StoreIntermParticles(0)
{
  std::cout << "voftopo instance created" << std::endl;
  this->SetNumberOfInputPorts(3);
  this->Controller = vtkMPIController::New();
  this->Boundaries = vtkPolyData::New();
  this->VofGrid[0] = vtkRectilinearGrid::New();
  this->VofGrid[1] = vtkRectilinearGrid::New();
  this->VelocityGrid[0] = vtkRectilinearGrid::New();
  this->VelocityGrid[1] = vtkRectilinearGrid::New();
  IntermParticles.clear();
  IntermParticlesTimeStamps.clear();
  IntermParticleIds.clear();
  IntermParticleProcs.clear();
  g_emf0 = EMF0;
  g_emf1 = EMF1;
}

//----------------------------------------------------------------------------
vtkVofTopo::~vtkVofTopo()
{
  if (Seeds != 0) {
    Seeds->Delete();
  }
  this->Controller->Delete();
  this->Boundaries->Delete();
  this->VofGrid[0]->Delete();
  this->VofGrid[1]->Delete();
  this->VelocityGrid[0]->Delete();
  this->VelocityGrid[1]->Delete();
  std::cout << "voftopo instance destroyed" << std::endl;
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
  if (port == 2) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  }  
  return 1;
}

//----------------------------------------------------------------------------
void vtkVofTopo::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void getFluidExtent(vtkRectilinearGrid *vof, const int globalExtent[6],
		    const int numGhostLevels, int extent[6])
{
  vof->GetExtent(extent);
  extent[0] = extent[0] + (extent[0] > globalExtent[0] ? numGhostLevels : 0);
  extent[1] = extent[1] - (extent[1] < globalExtent[1] ? numGhostLevels : 0);
  extent[2] = extent[2] + (extent[2] > globalExtent[2] ? numGhostLevels : 0);
  extent[3] = extent[3] - (extent[3] < globalExtent[3] ? numGhostLevels : 0);
  extent[4] = extent[4] + (extent[4] > globalExtent[4] ? numGhostLevels : 0);
  extent[5] = extent[5] - (extent[5] < globalExtent[5] ? numGhostLevels : 0);

  int index;
  vtkDataArray *vofArray = vof->GetCellData()->GetArray("Data", index);
  if (index == -1) {
    std::cout << "KURWA" << std::endl;
  }
  int nodeRes[3];
  vof->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

  int imin = extent[0];
  int imax = extent[1];
  int jmin = extent[2];
  int jmax = extent[3];
  int kmin = extent[4];
  int kmax = extent[5];

  extent[0] = extent[2] = extent[4] = 100000;
  extent[1] = extent[3] = extent[5] = -100000;
  
  for (int k = kmin; k < kmax; ++k) {
    for (int j = jmin; j < jmax; ++j) {
      for (int i = imin; i < imax; ++i) {

	int idx = i + j*cellRes[0] + k*cellRes[0]*cellRes[1];
	float f = vofArray->GetComponent(idx, 0);

	if (f > g_emf0) {
	  if (extent[0] > i) extent[0] = i;
	  if (extent[1] < i) extent[1] = i;
	  if (extent[2] > j) extent[2] = j;
	  if (extent[3] < j) extent[3] = j;
	  if (extent[4] > k) extent[4] = k;
	  if (extent[5] < k) extent[5] = k;
	}
      }
    }
  }    
}

void vtkVofTopo::InterpolateField(vtkRectilinearGrid *vof[2],
				  vtkRectilinearGrid *velocity[2],
				  vtkRectilinearGrid *intVof,
				  vtkRectilinearGrid *intVelocity,
				  const float a)
{
  float dt = InputTimeValues[TimestepT1] - InputTimeValues[TimestepT0];
  if (TimeStepDelta != 0.0) {
    dt = TimeStepDelta;
  }
  dt *= Incr;

  int nodeRes[3];
  vof[0]->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  int optExtent[6];
  {
    int fluidExtent0[6];
    getFluidExtent(vof[0], GlobalExtent, NumGhostLevels, fluidExtent0);

    int fluidExtent1[6];
    getFluidExtent(vof[1], GlobalExtent, NumGhostLevels, fluidExtent1);

    optExtent[0] = std::min(fluidExtent0[0],fluidExtent1[0]);
    optExtent[1] = std::max(fluidExtent0[1],fluidExtent1[1]);
    optExtent[2] = std::min(fluidExtent0[2],fluidExtent1[2]);
    optExtent[3] = std::max(fluidExtent0[3],fluidExtent1[3]);
    optExtent[4] = std::min(fluidExtent0[4],fluidExtent1[4]);
    optExtent[5] = std::max(fluidExtent0[5],fluidExtent1[5]);
      
    optExtent[0] -= 4;
    optExtent[2] -= 4;
    optExtent[4] -= 4;

    optExtent[1] += 4;
    optExtent[3] += 4;
    optExtent[5] += 4;
  }
  
  std::vector<float4> particles;
  particles.clear();
  generateSeedPointsInCellCenters(VelocityGrid[0], Refinement, particles, optExtent, GlobalExtent, NumGhostLevels);

  std::vector<float4> particlesForward = particles;
  float t = InputTimeValues[TimestepT0];

  advectParticles(velocity, particlesForward, t, 1.0f,
		  InputTimeValues[TimestepT0],
		  InputTimeValues[TimestepT0] + (InputTimeValues[TimestepT1]-InputTimeValues[TimestepT0])/2.0f,
		  RK4NumSteps);
  
  std::vector<float4> particlesBackward = particles;
  t = InputTimeValues[TimestepT1];
  
  advectParticles(velocity, particlesBackward, t, -1.0f,
		  InputTimeValues[TimestepT0] + (InputTimeValues[TimestepT1]-InputTimeValues[TimestepT0])/2.0f,
		  InputTimeValues[TimestepT1], RK4NumSteps);

  intVof->CopyStructure(vof[0]);
  intVelocity->CopyStructure(velocity[0]);

  int index0,index1;
  vtkDataArray *vofArray0 = vof[0]->GetCellData()->GetArray("Data", index0);
  vtkDataArray *vofArray1 = vof[1]->GetCellData()->GetArray("Data", index1);
  vtkSmartPointer<vtkFloatArray> vofArray = vtkSmartPointer<vtkFloatArray>::New();
  vofArray->SetName("Data");
  vofArray->SetNumberOfComponents(1);
  vofArray->SetNumberOfTuples(vofArray0->GetNumberOfTuples());
  vofArray->FillComponent(0,0.0f);
  vtkDataArray *velocityArray0 = velocity[0]->GetCellData()->GetArray("Data", index0);
  vtkDataArray *velocityArray1 = velocity[1]->GetCellData()->GetArray("Data", index1);
  vtkSmartPointer<vtkFloatArray> velocityArray = vtkSmartPointer<vtkFloatArray>::New();
  velocityArray->SetName("Data");
  velocityArray->SetNumberOfComponents(3);
  velocityArray->SetNumberOfTuples(velocityArray0->GetNumberOfTuples());
  velocityArray->FillComponent(0,0.0f);
  velocityArray->FillComponent(1,0.0f);
  velocityArray->FillComponent(2,0.0f);
  
  //#pragma omp parallel for
  for (int i = 0; i < particles.size(); ++i) {

    int ijk[3];
    double pcoords[3];    
    double x[3] = {particles[i].x, particles[i].y, particles[i].z};    
    vof[0]->ComputeStructuredCoordinates(x, ijk, pcoords);    
    int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
    
    float f0 = vofArray0->GetComponent(idx,0);
    float3 v0 = make_float3(velocityArray0->GetComponent(idx,0),
			    velocityArray0->GetComponent(idx,1),
			    velocityArray0->GetComponent(idx,2));
    float f1 = vofArray1->GetComponent(idx,0);
    float3 v1 = make_float3(velocityArray1->GetComponent(idx,0),
			    velocityArray1->GetComponent(idx,1),
			    velocityArray1->GetComponent(idx,2));
    
    double x0[3] = {particlesForward[i].x, particlesForward[i].y, particlesForward[i].z};
    vof[0]->ComputeStructuredCoordinates(x0, ijk, pcoords);
    {
      int lx = ijk[0];
      int ly = ijk[1];
      int lz = ijk[2];
      float bx = pcoords[0] - 0.5;
      float by = pcoords[1] - 0.5;
      float bz = pcoords[2] - 0.5;

      if (pcoords[0] < 0.5) {
	lx -= 1;
	bx = pcoords[0] + 0.5;
      }
      if (pcoords[1] < 0.5) {
	ly -= 1;
	by = pcoords[1] + 0.5;
      }
      if (pcoords[2] < 0.5) {
	lz -= 1;
	bz = pcoords[2] + 0.5;
      }
   
      int ux = lx+1;
      int uy = ly+1;
      int uz = lz+1;
    
      if (lx < 0) lx = 0;
      if (ly < 0) ly = 0;
      if (lz < 0) lz = 0;
      if (ux > cellRes[0]-1) ux = cellRes[0]-1;
      if (uy > cellRes[1]-1) uy = cellRes[1]-1;
      if (uz > cellRes[2]-1) uz = cellRes[2]-1;

      float contr[8] = {(1.0f-bx)*(1.0f-by)*(1.0f-bz),
			(bx)     *(1.0f-by)*(1.0f-bz),
			(1.0f-bx)*(by)     *(1.0f-bz),
			(bx)     *(by)     *(1.0f-bz),
			(1.0f-bx)*(1.0f-by)*(bz),
			(bx)     *(1.0f-by)*(bz),
			(1.0f-bx)*(by)     *(bz),
			(bx)     *(by)     *(bz)};

      int lzslab = lz*cellRes[0]*cellRes[1];
      int uzslab = uz*cellRes[0]*cellRes[1];
      int lyr = ly*cellRes[0];
      int uyr = uy*cellRes[0];
    
      int id[8] = {lx + lyr + lzslab,
		   ux + lyr + lzslab,
		   lx + uyr + lzslab,
		   ux + uyr + lzslab,
		   lx + lyr + uzslab,
		   ux + lyr + uzslab,
		   lx + uyr + uzslab,
		   ux + uyr + uzslab};

      for (int j = 0; j < 8; ++j) {
	vofArray->SetComponent(id[j],0,vofArray->GetComponent(id[j],0)+contr[j]*f0);
	velocityArray->SetComponent(id[j],0,velocityArray->GetComponent(id[j],0)+contr[j]*v0.x);
	velocityArray->SetComponent(id[j],1,velocityArray->GetComponent(id[j],1)+contr[j]*v0.y);
	velocityArray->SetComponent(id[j],2,velocityArray->GetComponent(id[j],2)+contr[j]*v0.z);
      }
    }
        
    double x1[3] = {particlesBackward[i].x, particlesBackward[i].y, particlesBackward[i].z};
    vof[0]->ComputeStructuredCoordinates(x1, ijk, pcoords);
    {
      int lx = ijk[0];
      int ly = ijk[1];
      int lz = ijk[2];
      float bx = pcoords[0] - 0.5;
      float by = pcoords[1] - 0.5;
      float bz = pcoords[2] - 0.5;

      if (pcoords[0] < 0.5) {
	lx -= 1;
	bx = pcoords[0] + 0.5;
      }
      if (pcoords[1] < 0.5) {
	ly -= 1;
	by = pcoords[1] + 0.5;
      }
      if (pcoords[2] < 0.5) {
	lz -= 1;
	bz = pcoords[2] + 0.5;
      }
   
      int ux = lx+1;
      int uy = ly+1;
      int uz = lz+1;
    
      if (lx < 0) lx = 0;
      if (ly < 0) ly = 0;
      if (lz < 0) lz = 0;
      if (ux > cellRes[0]-1) ux = cellRes[0]-1;
      if (uy > cellRes[1]-1) uy = cellRes[1]-1;
      if (uz > cellRes[2]-1) uz = cellRes[2]-1;

      float contr[8] = {(1.0f-bx)*(1.0f-by)*(1.0f-bz),
			(bx)     *(1.0f-by)*(1.0f-bz),
			(1.0f-bx)*(by)     *(1.0f-bz),
			(bx)     *(by)     *(1.0f-bz),
			(1.0f-bx)*(1.0f-by)*(bz),
			(bx)     *(1.0f-by)*(bz),
			(1.0f-bx)*(by)     *(bz),
			(bx)     *(by)     *(bz)};

      int lzslab = lz*cellRes[0]*cellRes[1];
      int uzslab = uz*cellRes[0]*cellRes[1];
      int lyr = ly*cellRes[0];
      int uyr = uy*cellRes[0];
    
      int id[8] = {lx + lyr + lzslab,
		   ux + lyr + lzslab,
		   lx + uyr + lzslab,
		   ux + uyr + lzslab,
		   lx + lyr + uzslab,
		   ux + lyr + uzslab,
		   lx + uyr + uzslab,
		   ux + uyr + uzslab};

      for (int j = 0; j < 8; ++j) {
	vofArray->SetComponent(id[j],0,vofArray->GetComponent(id[j],0)+contr[j]*f1);
	velocityArray->SetComponent(id[j],0,velocityArray->GetComponent(id[j],0)+contr[j]*v1.x);
	velocityArray->SetComponent(id[j],1,velocityArray->GetComponent(id[j],1)+contr[j]*v1.y);
	velocityArray->SetComponent(id[j],2,velocityArray->GetComponent(id[j],2)+contr[j]*v1.z);
      }
    }
  }
  const int numCells = cellRes[0]*cellRes[1]*cellRes[2];
  for (int idx = 0; idx < numCells; ++idx) {
    float f = vofArray->GetComponent(idx,0);
    vofArray->SetComponent(idx,0,f/2.0f);

    float3 v = make_float3(velocityArray->GetComponent(idx,0),
  			   velocityArray->GetComponent(idx,1),
  			   velocityArray->GetComponent(idx,2));
    v /= 2.0f;
    velocityArray->SetComponent(idx,0,v.x);
    velocityArray->SetComponent(idx,1,v.y);
    velocityArray->SetComponent(idx,2,v.z);
  }
  
  intVof->GetCellData()->AddArray(vofArray);
  intVelocity->GetCellData()->AddArray(velocityArray);
}

// void vtkVofTopo::InterpolateField(vtkRectilinearGrid *vof[2],
// 				  vtkRectilinearGrid *velocity[2],
// 				  vtkRectilinearGrid *intVof,
// 				  vtkRectilinearGrid *intVelocity,
// 				  const float a)
// {
//   int nodeRes[3];
//   vof[0]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

//   int optExtent[6];

//   {
//     int fluidExtent0[6];
//     getFluidExtent(vof[0], GlobalExtent, NumGhostLevels, fluidExtent0);

//     int fluidExtent1[6];
//     getFluidExtent(vof[1], GlobalExtent, NumGhostLevels, fluidExtent1);

//     optExtent[0] = std::min(fluidExtent0[0],fluidExtent1[0]);
//     optExtent[1] = std::max(fluidExtent0[1],fluidExtent1[1]);
//     optExtent[2] = std::min(fluidExtent0[2],fluidExtent1[2]);
//     optExtent[3] = std::max(fluidExtent0[3],fluidExtent1[3]);
//     optExtent[4] = std::min(fluidExtent0[4],fluidExtent1[4]);
//     optExtent[5] = std::max(fluidExtent0[5],fluidExtent1[5]);
      
//     optExtent[0] -= 4;
//     optExtent[2] -= 4;
//     optExtent[4] -= 4;

//     optExtent[1] += 4;
//     optExtent[3] += 4;
//     optExtent[5] += 4;
//   }
  
//   std::vector<float4> particles;
//   particles.clear();
//   generateSeedPointsInCellCenters(VelocityGrid[0], Refinement, particles, optExtent, GlobalExtent, NumGhostLevels);
  
  
//   float dt = InputTimeValues[TimestepT1] - InputTimeValues[TimestepT0];
//   if (TimeStepDelta != 0.0) {
//     dt = TimeStepDelta;
//   }
  
//   dt *= Incr;

//   std::vector<float4> particlesForward = particles;
//   float t = (1.0f-a)*InputTimeValues[TimestepT1] + a*InputTimeValues[TimestepT0];

//   advectParticles(velocity, particlesForward, t, 1.0f,
// 		  InputTimeValues[TimestepT0], InputTimeValues[TimestepT1], RK4NumSteps);

//   // {
//   //   vtkPolyData *outParticles = vtkPolyData::New();
//   //   vtkPoints *ppoints = vtkPoints::New();
//   //   ppoints->SetNumberOfPoints(particlesForward.size());
    
//   //   for (int i = 0; i < particlesForward.size(); ++i) {

//   //     float p[3] = {particlesForward[i].x, particlesForward[i].y, particlesForward[i].z};
//   //     ppoints->SetPoint(i, p);
//   //   }
//   //   outParticles->SetPoints(ppoints);

//   //   vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//   //   writer->SetInputData(outParticles);
//   //   writer->SetFileName("/tmp/intPartForward.vtp");
//   //   writer->Write();
//   // }
 
//   std::vector<float4> particlesBackward = particles;
//   advectParticles(velocity, particlesBackward, t, -1.0f,
// 		  InputTimeValues[TimestepT0], InputTimeValues[TimestepT1], RK4NumSteps);

//   // {
//   //   vtkPolyData *outParticles = vtkPolyData::New();
//   //   vtkPoints *ppoints = vtkPoints::New();
//   //   ppoints->SetNumberOfPoints(particlesBackward.size());
    
//   //   for (int i = 0; i < particlesBackward.size(); ++i) {

//   //     float p[3] = {particlesBackward[i].x, particlesBackward[i].y, particlesBackward[i].z};
//   //     ppoints->SetPoint(i, p);
//   //   }
//   //   outParticles->SetPoints(ppoints);

//   //   vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//   //   writer->SetInputData(outParticles);
//   //   writer->SetFileName("/tmp/intPartBackward.vtp");
//   //   writer->Write();
//   // }

//   intVof->CopyStructure(vof[0]);
//   intVelocity->CopyStructure(velocity[0]);

//   int index0,index1;
//   vtkDataArray *vofArray0 = vof[0]->GetCellData()->GetArray("Data", index0);
//   vtkDataArray *vofArray1 = vof[1]->GetCellData()->GetArray("Data", index1);
//   vtkSmartPointer<vtkFloatArray> vofArray = vtkSmartPointer<vtkFloatArray>::New();
//   vofArray->SetName("Data");
//   vofArray->SetNumberOfComponents(1);
//   vofArray->SetNumberOfTuples(vofArray0->GetNumberOfTuples());
//   vofArray->FillComponent(0,0.0f);
//   vtkDataArray *velocityArray0 = velocity[0]->GetCellData()->GetArray("Data", index0);
//   vtkDataArray *velocityArray1 = velocity[1]->GetCellData()->GetArray("Data", index1);
//   vtkSmartPointer<vtkFloatArray> velocityArray = vtkSmartPointer<vtkFloatArray>::New();
//   velocityArray->SetName("Data");
//   velocityArray->SetNumberOfComponents(3);
//   velocityArray->SetNumberOfTuples(velocityArray0->GetNumberOfTuples());
//   velocityArray->FillComponent(0,0.0f);
//   velocityArray->FillComponent(1,0.0f);
//   velocityArray->FillComponent(2,0.0f);
  
// #pragma omp parallel for
//   for (int i = 0; i < particles.size(); ++i) {

//     int ijk[3];
//     double pcoords[3];
    
//     double x0[3] = {particlesBackward[i].x, particlesBackward[i].y, particlesBackward[i].z};    
//     vof[0]->ComputeStructuredCoordinates(x0, ijk, pcoords);
//     float f0 = interpolateScaCellBasedData(vofArray0, cellRes, ijk, pcoords);
//     float3 v0 = interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords);
    
//     double x1[3] = {particlesForward[i].x, particlesForward[i].y, particlesForward[i].z};
//     vof[1]->ComputeStructuredCoordinates(x1, ijk, pcoords);
//     float f1 = interpolateScaCellBasedData(vofArray1, cellRes, ijk, pcoords);
//     float3 v1 = interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords);
    
//     double x[3] = {particles[i].x, particles[i].y, particles[i].z};
//     vof[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
//     float f = (1.0f-a)*f0 + a*f1;
//     float3 v = (1.0f-a)*v0 + a*v1;

//     int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];    
//     vofArray->SetValue(idx,f);
//     velocityArray->SetTuple3(idx,v.x,v.y,v.z);
//   }
//   intVof->GetCellData()->AddArray(vofArray);
//   intVelocity->GetCellData()->AddArray(velocityArray);

//   // {
//   //   vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
//   //   writer->SetInputData(intVof);
//   //   writer->SetFileName("/tmp/intVof.vtr");
//   //   writer->Write();
//   // }
//   // {
//   //   vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
//   //   writer->SetInputData(intVelocity);
//   //   writer->SetFileName("/tmp/intVelocity.vtr");
//   //   writer->Write();
//   // }
// }

void vtkVofTopo::TransferIntermParticlesToSeeds(std::vector<std::vector<float4>> &particles,
						std::vector<std::vector<int>> &ids,
						std::vector<std::vector<short>> &procs)
{
  if (Controller->GetCommunicator() == 0) {
    return;
  }
  const int numProcesses = Controller->GetNumberOfProcesses();
  const int processId = Controller->GetLocalProcessId();
  const int numTimeSteps = particles.size();
  
  for (int i = 0; i < numTimeSteps; ++i) {

    // one vector for each side of the process
    std::vector<std::vector<float4>> particlesToSend(numProcesses);
    std::vector<std::vector<int> > idsToSend(numProcesses);
  
    for (int j = 0; j < numProcesses; ++j) {
      particlesToSend[j].resize(0);
      idsToSend[j].resize(0);
    }

    std::vector<float4> particlesToKeep(0);
    std::vector<int> idsToKeep(0);

    for (int j = 0; j < particles[i].size(); ++j) {

      int procId = procs[i][j];
      if (processId != procId) {
	particlesToSend[procId].push_back(particles[i][j]);
	idsToSend[procId].push_back(ids[i][j]);
      }
      else {
	particlesToKeep.push_back(particles[i][j]);
	idsToKeep.push_back(ids[i][j]);
      }
    }
 
    particles[i] = particlesToKeep;
    ids[i] = idsToKeep;

    std::vector<float4> particlesToRecv;
    std::vector<int> idsToRecv;
    sendData(particlesToSend, particlesToRecv, numProcesses, Controller);
    sendData(idsToSend, idsToRecv, numProcesses, Controller);

    // insert the paricles that are within the domain
    for (int j = 0; j < particlesToRecv.size(); ++j) {
      particles[i].push_back(particlesToRecv[j]);
      ids[i].push_back(idsToRecv[j]);
    }

    Controller->GetCommunicator()->Barrier();    
  }
}

