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
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <unordered_set>
#include <array>
#include <string>
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkXMLPolyDataWriter.h"


vtkStandardNewMacro(vtkVofTopo);


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

// void writeData(vtkPolyData *data, const int blockId,
// 	       const int processId, const std::string path)
// {
//   vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//   writer->SetInputData(data);
//   std::string outFile = path;
//   outFile += std::to_string(blockId) + std::string("_") + std::to_string(processId);
//   outFile += ".vtp";
//   writer->SetFileName(outFile.c_str());
//   writer->SetDataModeToBinary();
//   writer->Write();
// }

//----------------------------------------------------------------------------
int vtkVofTopo::RequestData(vtkInformation *request,
			    vtkInformationVector **inputVector,
			    vtkInformationVector *outputVector)
{
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

  if (TimestepT0 == TimestepT1 && Controller->GetCommunicator() != 0) {
    // find neighbor processes and global domain bounds
    GetGlobalContext(inInfoVof);
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

      InitVelocities(VelocityGrid[0]);

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
      AdvectParticles(VofGrid, VelocityGrid);
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
	TransferLabelsToSeeds(particleLabels);

	// Transfer seed points from neighbors ---------------------------------
	vtkPolyData *boundarySeeds = vtkPolyData::New();
	if (Controller->GetCommunicator() != 0) {
	  ExchangeBoundarySeedPoints(boundarySeeds);
	}
	
	// Stage VI ------------------------------------------------------------
	GenerateBoundaries(Boundaries, boundarySeeds);

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
	output->SetBlock(1, particles);

	output->SetBlock(2, Boundaries);
	output->SetBlock(0, Seeds);

	// writeData(Seeds, 0, Controller->GetLocalProcessId(), "out2/out2_");
	// writeData(particles, 1, Controller->GetLocalProcessId(), "out2/out2_");
	// writeData(Boundaries, 2, Controller->GetLocalProcessId(), "out2/out2_");
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
  //  vtkSmartPointer<vtkPoints> seedPoints = vtkSmartPointer<vtkPoints>::New();
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
void vtkVofTopo::InitVelocities(vtkRectilinearGrid *velocity)
{
  Velocities.clear();
  vtkPoints *seedPoints = Seeds->GetPoints();
  Velocities.resize(seedPoints->GetNumberOfPoints());

  initVelocities(velocity, Particles, Velocities);
}
void computeBoundsWithoutGhost(double globalBounds[6], double localBounds[6], 
			       vtkDataArray *xCoordinates, 
			       vtkDataArray *yCoordinates, 
			       vtkDataArray *zCoordinates,
			       int numGhostLevels,
			       double boundsNoGhosts[6])
{
  boundsNoGhosts[0] = localBounds[0];
  boundsNoGhosts[1] = localBounds[1];
  boundsNoGhosts[2] = localBounds[2];
  boundsNoGhosts[3] = localBounds[3];
  boundsNoGhosts[4] = localBounds[4];
  boundsNoGhosts[5] = localBounds[5];

  if (localBounds[0] != globalBounds[0]) {
    boundsNoGhosts[0] = xCoordinates->GetComponent(numGhostLevels,0);
  }
  if (localBounds[1] != globalBounds[1]) {
    int numCoords = xCoordinates->GetNumberOfTuples();
    boundsNoGhosts[1] = xCoordinates->GetComponent(numCoords-1-numGhostLevels,0);
  }
  if (localBounds[2] != globalBounds[2]) {
    boundsNoGhosts[2] = yCoordinates->GetComponent(numGhostLevels,0);
  }
  if (localBounds[3] != globalBounds[3]) {
    int numCoords = yCoordinates->GetNumberOfTuples();
    boundsNoGhosts[3] = yCoordinates->GetComponent(numCoords-1-numGhostLevels,0);
  }
  if (localBounds[4] != globalBounds[4]) {
    boundsNoGhosts[4] = zCoordinates->GetComponent(numGhostLevels,0);
  }
  if (localBounds[5] != globalBounds[5]) {
    int numCoords = zCoordinates->GetNumberOfTuples();
    boundsNoGhosts[5] = zCoordinates->GetComponent(numCoords-1-numGhostLevels,0);
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

  findGlobalExtent(AllExtents, GlobalExtent);

  // reduce extent to one without ghost cells
  if (LocalExtents[0] > GlobalExtent[0])
    LocalExtents[0] += NumGhostLevels;
  if (LocalExtents[1] < GlobalExtent[1])
    LocalExtents[1] -= NumGhostLevels;
  if (LocalExtents[2] > GlobalExtent[2])
    LocalExtents[2] += NumGhostLevels;
  if (LocalExtents[3] < GlobalExtent[3])
    LocalExtents[3] -= NumGhostLevels;
  if (LocalExtents[4] > GlobalExtent[4])
    LocalExtents[4] += NumGhostLevels;
  if (LocalExtents[5] < GlobalExtent[5])
    LocalExtents[5] -= NumGhostLevels;

  // send extents again
  Controller->AllGatherV(&LocalExtents[0], &AllExtents[0], NUM_SIDES, &RecvLengths[0], &RecvOffsets[0]);

  NeighborProcesses.clear();
  NeighborProcesses.resize(NUM_SIDES);
  findNeighbors(LocalExtents, GlobalExtent, AllExtents, NeighborProcesses);

  NumNeighbors = 0;
  for (int i = 0; i < NeighborProcesses.size(); ++i) {
    NumNeighbors += NeighborProcesses[i].size();
  }

  // find domain bounds
  inputVof->GetBounds(&LocalBounds[0]);
  std::vector<double> AllBounds(NUM_SIDES*numProcesses);
  Controller->AllGatherV(&LocalBounds[0], &AllBounds[0], 6, &RecvLengths[0], &RecvOffsets[0]);
  findGlobalBounds(AllBounds, GlobalBounds);

  computeBoundsWithoutGhost(GlobalBounds, LocalBounds, 
			    inputVof->GetXCoordinates(), 
			    inputVof->GetYCoordinates(), 
			    inputVof->GetZCoordinates(),
			    NumGhostLevels,
			    BoundsNoGhosts);
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
  
  advectParticles(vof[1], velocity[1], Particles, Velocities, dt, Uncertainty);
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
  std::vector<std::vector<float4>> particlesToSend(numProcesses);
  std::vector<std::vector<float4>> velocitiesToSend(numProcesses);
  std::vector<std::vector<int> > particleIdsToSend(numProcesses);
  std::vector<std::vector<short> > particleProcsToSend(numProcesses);
  std::vector<std::vector<float>> uncertaintyToSend(numProcesses);
  
  for (int i = 0; i < numProcesses; ++i) {
    particlesToSend[i].resize(0);
    velocitiesToSend[i].resize(0);
    particleIdsToSend[i].resize(0);
    particleProcsToSend[i].resize(0);
    uncertaintyToSend[i].resize(0);
  }

  std::vector<float4>::iterator it;
  std::vector<float4> particlesToKeep;
  std::vector<float4> velocitiesToKeep;
  std::vector<int> particleIdsToKeep;
  std::vector<short> particleProcsToKeep;
  std::vector<float> uncertaintyToKeep;

  for (int i = 0; i < Particles.size(); ++i) {

    int bound = outOfBounds(Particles[i], BoundsNoGhosts, GlobalBounds);
    if (bound > -1) {
      for (int j = 0; j < NeighborProcesses[bound].size(); ++j) {

  	int neighborId = NeighborProcesses[bound][j];
	if (neighborId != processId) {
	  particlesToSend[neighborId].push_back(Particles[i]);
	  velocitiesToSend[neighborId].push_back(Velocities[i]);
	  particleIdsToSend[neighborId].push_back(ParticleIds[i]);
	  particleProcsToSend[neighborId].push_back(ParticleProcs[i]);
	  uncertaintyToSend[neighborId].push_back(Uncertainty[i]);
	}
      }
    }
    else {
      particlesToKeep.push_back(Particles[i]);
      velocitiesToKeep.push_back(Velocities[i]);
      particleIdsToKeep.push_back(ParticleIds[i]);
      particleProcsToKeep.push_back(ParticleProcs[i]);
      uncertaintyToKeep.push_back(Uncertainty[i]);
    }
  }
 
  Particles = particlesToKeep;
  Velocities = velocitiesToKeep;
  ParticleIds = particleIdsToKeep;
  ParticleProcs = particleProcsToKeep;
  Uncertainty = uncertaintyToKeep;

  std::vector<float4> particlesToRecv;
  std::vector<float4> velocitiesToRecv;
  std::vector<int> particleIdsToRecv;
  std::vector<short> particleProcsToRecv;
  std::vector<float> uncertaintyToRecv;
  sendData(particlesToSend, particlesToRecv, numProcesses, Controller);
  sendData(velocitiesToSend, velocitiesToRecv, numProcesses, Controller);
  sendData(particleIdsToSend, particleIdsToRecv, numProcesses, Controller);
  sendData(particleProcsToSend, particleProcsToRecv, numProcesses, Controller);
  sendData(uncertaintyToSend, uncertaintyToRecv, numProcesses, Controller);

  // insert the paricles that are within the domain
  for (int i = 0; i < particlesToRecv.size(); ++i) {
    int within = withinBounds(particlesToRecv[i], BoundsNoGhosts);
    if (within) {
      Particles.push_back(particlesToRecv[i]);
      Velocities.push_back(velocitiesToRecv[i]);
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

    std::vector<std::vector<float>> labelsToSend(numProcesses);
    std::vector<std::vector<int>> idsToSend(numProcesses);
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
void vtkVofTopo::GenerateBoundaries(vtkPolyData *boundaries, vtkPolyData *boundarySeeds)
{
  vtkPoints *points = Seeds->GetPoints();
  vtkFloatArray *labels = vtkFloatArray::
    SafeDownCast(Seeds->GetPointData()->GetArray("Labels"));

  vtkPoints *boundarySeedPoints = boundarySeeds->GetPoints();
  vtkFloatArray *boundarySeedLabels = vtkFloatArray::
    SafeDownCast(boundarySeeds->GetPointData()->GetArray("Labels"));
  
  generateBoundaries(points, boundarySeedPoints,
		     labels, boundarySeedLabels,
		     this->VofGrid[1], boundaries, this->Refinement);
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
  int extent[6];
  this->VofGrid[1]->GetExtent(extent);

  int nodeRes[3] = {this->VofGrid[1]->GetXCoordinates()->GetNumberOfTuples(), 
		    this->VofGrid[1]->GetYCoordinates()->GetNumberOfTuples(), 
		    this->VofGrid[1]->GetZCoordinates()->GetNumberOfTuples()};
  int cellRes[3] = {nodeRes[0]-1,nodeRes[1]-1,nodeRes[2]-1};
  
  int imin = extent[0] > GlobalExtent[0] ? NumGhostLevels : 0;
  int imax = extent[1] < GlobalExtent[1] ? cellRes[0]-NumGhostLevels : cellRes[0];
  int jmin = extent[2] > GlobalExtent[2] ? NumGhostLevels : 0;
  int jmax = extent[3] < GlobalExtent[3] ? cellRes[1]-NumGhostLevels : cellRes[1];
  int kmin = extent[4] > GlobalExtent[4] ? NumGhostLevels : 0;
  int kmax = extent[5] < GlobalExtent[5] ? cellRes[2]-NumGhostLevels : cellRes[2];

  int bids[6] = {extent[0] > GlobalExtent[0] ? NumGhostLevels : -1,
		 extent[1] < GlobalExtent[1] ? cellRes[0]-NumGhostLevels-1 : -1,
		 extent[2] > GlobalExtent[2] ? NumGhostLevels : -1,
		 extent[3] < GlobalExtent[3] ? cellRes[1]-NumGhostLevels-1 : -1,
		 extent[4] > GlobalExtent[4] ? NumGhostLevels : -1,
		 extent[5] < GlobalExtent[5] ? cellRes[2]-NumGhostLevels-1 : -1};

  std::set<std::array<int,3>> boundaryCells;
  boundaryCells.clear();

  if (bids[0] > -1) { // left
    const int i = bids[0];
    for (int k = kmin; k < kmax; ++k) {
      for (int j = jmin; j < jmax; ++j) {
	boundaryCells.emplace(std::array<int,3>{i,j,k});
      }
    }
  }
  if (bids[1] > -1) { // right
    const int i = bids[1];
    for (int k = kmin; k < kmax; ++k) {
      for (int j = jmin; j < jmax; ++j) {
	boundaryCells.emplace(std::array<int,3>{i,j,k});
      }
    }
  }
  if (bids[2] > -1) { // bottom
    const int j = bids[2];
    for (int k = kmin; k < kmax; ++k) {
      for (int i = imin; i < imax; ++i) {
	boundaryCells.emplace(std::array<int,3>{i,j,k});
      }
    }
  }
  if (bids[3] > -1) { // top
    const int j = bids[3];
    for (int k = kmin; k < kmax; ++k) {
      for (int i = imin; i < imax; ++i) {
	boundaryCells.emplace(std::array<int,3>{i,j,k});
      }
    }
  }
  if (bids[4] > -1) { // back
    const int k = bids[4];
    for (int j = jmin; j < jmax; ++j) {
      for (int i = imin; i < imax; ++i) {
	boundaryCells.emplace(std::array<int,3>{i,j,k});
      }
    }
  }
  if (bids[5] > -1) { // front
    const int k = bids[5];
    for (int j = jmin; j < jmax; ++j) {
      for (int i = imin; i < imax; ++i) {
	boundaryCells.emplace(std::array<int,3>{i,j,k});
      }
    }
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
    int insideGrid = isInside(ijk,cellRes,NumGhostLevels,1);
    if (inside && insideGrid) {
    
      boundarySeedPoints->InsertNextPoint(x);
      boundarySeedLabels->InsertNextTuple1(labelsToRecv[i]);
    }
  }

  // if (pointsToRecv.size() > 0) {
  //   std::cout << "be: "
  // 	      << be[0] << "-" << be[1] << " x "
  // 	      << be[2] << "-" << be[3] << " x "
  // 	      << be[4] << "-" << be[5] << std::endl
  // 	      << "extent: "
  // 	      << extent[0] << "-" << extent[1] << " x "
  // 	      << extent[2] << "-" << extent[3] << " x "
  // 	      << extent[4] << "-" << extent[5] << std::endl;
  // }
  
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
  SeedPointsProvided(false)
{
  this->SetNumberOfInputPorts(3);
  this->Controller = vtkMPIController::New();
  this->Boundaries = vtkPolyData::New();
  this->VofGrid[0] = vtkRectilinearGrid::New();
  this->VofGrid[1] = vtkRectilinearGrid::New();
  this->VelocityGrid[0] = vtkRectilinearGrid::New();
  this->VelocityGrid[1] = vtkRectilinearGrid::New();
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
