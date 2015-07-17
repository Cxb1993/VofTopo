#include "vtkVofTopo.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
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
  CurrentTimeStep(0)
{
  this->SetNumberOfInputPorts(2);
}

//----------------------------------------------------------------------------
vtkVofTopo::~vtkVofTopo()
{
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
      vtkWarningMacro(<<"InitTimeStep out of range; setting to " << InitTimeStep);
    }
    if (InitTimeStep > InputTimeValues.size()-1) {
      InitTimeStep = InputTimeValues.size()-1;
      vtkWarningMacro(<<"InitTimeStep out of range; setting to " << InitTimeStep);
    }
    if (TargetTimeStep > InputTimeValues.size()-1) {
      TargetTimeStep = InputTimeValues.size()-1;
      vtkWarningMacro(<<"TargetTimeStep out of range; setting to " << TargetTimeStep);
    }
    if (TargetTimeStep < InitTimeStep) {
      vtkErrorMacro(<<"TargetTimeStep smaller than InitTimeStep");
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
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  // inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  // outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);
  
  if(FirstIteration) {

    if (IterType == IterateTarget) {

      double targetTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      if(targetTime > InputTimeValues.back()) {
      	targetTime = InputTimeValues.back();
      }
      TargetTimeStep = findClosestTimeStep(targetTime, InputTimeValues);
      if (TargetTimeStep < 0) {
	TargetTimeStep = 0;
	vtkWarningMacro(<<"TargetTimeStep out of range; setting to " << TargetTimeStep);
      }
      if (TargetTimeStep > InputTimeValues.size()-1) {
	TargetTimeStep = InputTimeValues.size()-1;
	vtkWarningMacro(<<"TargetTimeStep out of range; setting to " << TargetTimeStep);
      }
      CurrentTimeStep = InitTimeStep;
    }

  }
  if (InitTimeStep <= TargetTimeStep) {
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

  return 1;
}

//----------------------------------------------------------------------------
void vtkVofTopo::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
