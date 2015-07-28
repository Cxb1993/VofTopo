#include "vtkVofTopo.h"
#include "vofTopology.h"

#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkMPIController.h"
#include "vtkMPICommunicator.h"
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
      float dt = InputTimeValues[CurrentTimeStep+1] - InputTimeValues[CurrentTimeStep];
      advectParticles(inputVof, inputVelocity, Particles, dt);
      LastComputedTimeStep = CurrentTimeStep;
    }

    if (Controller->GetCommunicator() != 0) {
      ExchangeParticles();
    }

    bool finished = (CurrentTimeStep + 1) >= TargetTimeStep;
    // Stage III -------------------------------------------------------------
    if (finished) {
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      FirstIteration = true;

      vtkInformation *outInfo = outputVector->GetInformationObject(0);
      vtkPolyData *output =
      	vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
      GenerateOutputGeometry(output);
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
void vtkVofTopo::GenerateSeeds(vtkRectilinearGrid *inputVof)
{
  vtkSmartPointer<vtkPoints> seedPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> seedConnectivity = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkShortArray> seedCoords = vtkSmartPointer<vtkShortArray>::New();
  generateSeedPoints(inputVof, Refinement, seedPoints, seedConnectivity, seedCoords);
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
void vtkVofTopo::ExchangeParticles()
{
  int numProcesses = Controller->GetNumberOfProcesses();
  int processId = Controller->GetLocalProcessId();
  // one vector for each side of the process
  std::vector<std::vector<float4> > particlesToSend(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    particlesToSend[i].resize(0);
  }

  //------------------------------------------------------------------------
  std::vector<float4>::iterator it;
  std::vector<float4> particlesToKeep;

  int pidx = 0;
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

    ++pidx;
  }
  Particles = particlesToKeep;
  //------------------------------------------------------------------------

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

