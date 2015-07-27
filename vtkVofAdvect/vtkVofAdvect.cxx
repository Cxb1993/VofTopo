// TODO
// clean up code
// reuse particles from previous time step

#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkShortArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkRectilinearGrid.h"
#include "vtkPolyData.h"
#include "vtkMPIController.h"

#include "vtkVofAdvect.h"
#include "voftopo.h"

#include <iostream>
#include <vector>
#include <limits>
#include <set>
#include <map>

vtkStandardNewMacro(vtkVofAdvect);

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

  // if the bounds is on the global domain boundary, then don't mark the particle for sending
  int outOfBounds(const f3u1_t particle, const double bounds[6], const double globalBounds[6])
  {
    if (particle.x < bounds[0] && bounds[0] != globalBounds[0]) return 0;
    if (particle.x > bounds[1] && bounds[1] != globalBounds[1]) return 1;
    if (particle.y < bounds[2] && bounds[2] != globalBounds[2]) return 2;
    if (particle.y > bounds[3] && bounds[3] != globalBounds[3]) return 3;
    if (particle.z < bounds[4] && bounds[4] != globalBounds[4]) return 4;
    if (particle.z > bounds[5] && bounds[5] != globalBounds[5]) return 5;

    return -1;
  }

  int withinBounds(const f3u1_t particle, const double bounds[6])
  {
    if (particle.x < bounds[0]) return 0;
    if (particle.x > bounds[1]) return 0;
    if (particle.y < bounds[2]) return 0;
    if (particle.y > bounds[3]) return 0;
    if (particle.z < bounds[4]) return 0;
    if (particle.z > bounds[5]) return 0;

    return 1;
  }

  void findGlobalExtents(std::vector<int> &allExtents, 
			 int globalExtents[6])
  {
    globalExtents[0] = globalExtents[2] = globalExtents[4] = std::numeric_limits<int>::max();
    globalExtents[1] = globalExtents[3] = globalExtents[5] = - globalExtents[0];

    for (int i = 0; i < allExtents.size()/6; ++i) {
      if (globalExtents[0] > allExtents[i*6+0]) globalExtents[0] = allExtents[i*6+0];
      if (globalExtents[1] < allExtents[i*6+1]) globalExtents[1] = allExtents[i*6+1];
      if (globalExtents[2] > allExtents[i*6+2]) globalExtents[2] = allExtents[i*6+2];
      if (globalExtents[3] < allExtents[i*6+3]) globalExtents[3] = allExtents[i*6+3];
      if (globalExtents[4] > allExtents[i*6+4]) globalExtents[4] = allExtents[i*6+4];
      if (globalExtents[5] < allExtents[i*6+5]) globalExtents[5] = allExtents[i*6+5];
    }
  }

  void findGlobalBounds(std::vector<double> &allBounds, 
			double globalBounds[6])
  {
    globalBounds[0] = globalBounds[2] = globalBounds[4] = std::numeric_limits<double>::max();
    globalBounds[1] = globalBounds[3] = globalBounds[5] = - globalBounds[0];

    for (int i = 0; i < allBounds.size()/6; ++i) {
      if (globalBounds[0] > allBounds[i*6+0]) globalBounds[0] = allBounds[i*6+0];
      if (globalBounds[1] < allBounds[i*6+1]) globalBounds[1] = allBounds[i*6+1];
      if (globalBounds[2] > allBounds[i*6+2]) globalBounds[2] = allBounds[i*6+2];
      if (globalBounds[3] < allBounds[i*6+3]) globalBounds[3] = allBounds[i*6+3];
      if (globalBounds[4] > allBounds[i*6+4]) globalBounds[4] = allBounds[i*6+4];
      if (globalBounds[5] < allBounds[i*6+5]) globalBounds[5] = allBounds[i*6+5];
    }
  }

  void findNeighbors(const int myExtents[6], 
		     const int globalExtents[6], 
		     const std::vector<int> &allExtents,
		     std::vector<std::vector<int> > &neighbors)
  {
    const int numDims = 3;
    const int numSides = 6;

    for (int i = 0; i < numDims; ++i) {

      if (myExtents[i*2+0] > globalExtents[i*2+0]) { 
	for (int j = 0; j < allExtents.size()/numSides; ++j) {

	  if (myExtents[i*2+0] <= allExtents[j*numSides+i*2+1] &&
	      myExtents[i*2+1] > allExtents[j*numSides+i*2+1] &&
	      myExtents[((i+1)%3)*2+0] < allExtents[j*numSides+((i+1)%3)*2+1] &&
	      myExtents[((i+1)%3)*2+1] > allExtents[j*numSides+((i+1)%3)*2+0] &&
	      myExtents[((i+2)%3)*2+0] < allExtents[j*numSides+((i+2)%3)*2+1] &&
	      myExtents[((i+2)%3)*2+1] > allExtents[j*numSides+((i+2)%3)*2+0]) {

	    neighbors[i*2+0].push_back(j);
	  }
	}
      }
      if (myExtents[i*2+1] < globalExtents[i*2+1]) { 
	for (int j = 0; j < allExtents.size()/numSides; ++j) {

	  if (myExtents[i*2+1] >= allExtents[j*numSides+i*2+0] &&
	      myExtents[i*2+0] < allExtents[j*numSides+i*2+0] &&
	      myExtents[((i+1)%3)*2+0] < allExtents[j*numSides+((i+1)%3)*2+1] &&
	      myExtents[((i+1)%3)*2+1] > allExtents[j*numSides+((i+1)%3)*2+0] &&
	      myExtents[((i+2)%3)*2+0] < allExtents[j*numSides+((i+2)%3)*2+1] &&
	      myExtents[((i+2)%3)*2+1] > allExtents[j*numSides+((i+2)%3)*2+0]) {

	    neighbors[i*2+1].push_back(j);
	  }
	}
      }
    }
  }

  void getGridCenterCoords(vtkRectilinearGrid *grid,
			   vtkDataSetAttributes::AttributeTypes attrType,
			   vtkDataArray *coordCenters[3])
  {
    vtkDataArray *pointData = grid->GetPointData()->GetAttribute(attrType);
    vtkDataArray *cellData = grid->GetCellData()->GetAttribute(attrType);

    if (pointData == 0 && cellData != 0) {  // cell data
      vtkDataArray *coordNodes[3];
      coordNodes[0] = grid->GetXCoordinates();
      coordNodes[1] = grid->GetYCoordinates();
      coordNodes[2] = grid->GetZCoordinates();

      for (int c = 0; c < 3; ++c) {

	if (grid->GetXCoordinates()->GetDataType() == VTK_FLOAT) {
	  coordCenters[c] = vtkFloatArray::New();
	} else {
	  coordCenters[c] = vtkDoubleArray::New();
	}
	coordCenters[c]->SetNumberOfComponents(1);
	coordCenters[c]->SetNumberOfTuples(coordNodes[c]->GetNumberOfTuples()-1);
	for (int i = 0; i < coordCenters[c]->GetNumberOfTuples(); ++i) {
	  coordCenters[c]->SetComponent(i,0,(coordNodes[c]->GetComponent(0,i) + 
					     coordNodes[c]->GetComponent(0,i+1))/2.0f);
	}
      }
    }
    else if (pointData != 0 && cellData == 0) {  // point data
      coordCenters[0] = grid->GetXCoordinates();
      coordCenters[1] = grid->GetYCoordinates();
      coordCenters[2] = grid->GetZCoordinates();
    }
    else {
      std::cerr << "Ambiguity - both point data and cell data available" << std::endl;
    }
  }

  void getGridData(vtkRectilinearGrid *grid,
		   vtkDataSetAttributes::AttributeTypes attrType,
		   vtkDataArray **data)
  {
    vtkDataArray *pointData = grid->GetPointData()->GetAttribute(attrType);
    vtkDataArray *cellData = grid->GetCellData()->GetAttribute(attrType);

    if (pointData == 0 && cellData != 0) {  // cell data
      *data = cellData;
    }
    else if (pointData != 0 && cellData == 0) {  // point data
      *data = pointData;
    }
  }
  
  template<typename T> 
  void sendData(const std::vector<std::vector<T> > &dataToSend, std::vector<T> &dataToRecv,
		const int numProcesses, vtkMPIController *controller)
  {
    std::vector<T> allDataToSend(0);

    // flatten vector<vectr<T>> to vector<T>
    for (int i = 0; i < numProcesses; ++i) {
      for (int j = 0; j < dataToSend[i].size(); ++j) {
	allDataToSend.push_back(dataToSend[i][j]);
      }
    }

    std::vector<int> numDataToSend(numProcesses,0);
    for (int i = 0; i < numProcesses; ++i) {
      numDataToSend[i] = dataToSend[i].size();
    }

    // determine how many elements will be received from each process
    std::vector<int> numDataToRecv(numProcesses,0);
    std::vector<int> RecvLengths(numProcesses);
    std::vector<int> RecvOffsets(numProcesses);
    int numAllDataToRecv = 0;
    for (int i = 0; i < numProcesses; ++i) {

      controller->Scatter((int*)&numDataToSend[0], (int*)&numDataToRecv[i], 1, i);
      RecvOffsets[i] = numAllDataToRecv;
      RecvLengths[i] = numDataToRecv[i]*sizeof(T);
      numAllDataToRecv += numDataToRecv[i];
    }

    // send actual data
    dataToRecv.resize(numAllDataToRecv);

    std::vector<vtkIdType> SendLengths(numProcesses);
    std::vector<vtkIdType> SendOffsets(numProcesses);
    int offset = 0;
    for (int i = 0; i < numProcesses; ++i) {

      SendLengths[i] = numDataToSend[i]*sizeof(T);
      SendOffsets[i] = offset;
      offset += numDataToSend[i]*sizeof(T);
    }

    for (int i = 0; i < numProcesses; ++i) {
      controller->ScatterV((char*)&allDataToSend[0], 
			   (char*)&dataToRecv[RecvOffsets[i]], 
			   &SendLengths[0], &SendOffsets[0], RecvLengths[i], i);
    }
  }

};

//-----------------------------------------------------------------------------
vtkVofAdvect::vtkVofAdvect() :
  StartTime(0.0),
  TerminationTime(0.0),
  StartTimeStep(0),
  TerminationTimeStep(0),
  CurrentTimeStep(0),
  FirstIteration(true)
{
  this->SetNumberOfInputPorts(3);

  // get info on bounds from other processes
  this->Controller = vtkMPIController::New();
  int numProcesses = this->Controller->GetNumberOfProcesses();
}

//-----------------------------------------------------------------------------
vtkVofAdvect::~vtkVofAdvect()
{
  this->Controller->Delete();
}

//----------------------------------------------------------------------------
void vtkVofAdvect::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}

//----------------------------------------------------------------------------
void vtkVofAdvect::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
int vtkVofAdvect::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
  }
  if (port == 1) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
  }
  if (port == 2) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int vtkVofAdvect::RequestInformation(vtkInformation *vtkNotUsed(request),
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

    //clamp the default start time to be within the data time range
    if(this->StartTime < this->InputTimeValues[0]) {
      this->StartTime = this->InputTimeValues[0];
    }
    else if (this->StartTime > this->InputTimeValues.back()) {
      this->StartTime = this->InputTimeValues.back();
    }

    // check if user input is within time step range
    if (StartTimeStep < 0) {
      StartTimeStep = 0;
      vtkWarningMacro(<<"StartTimeStep out of range; setting to " << StartTimeStep);
    }
    if (StartTimeStep > InputTimeValues.size()-1) {
      StartTimeStep = InputTimeValues.size()-1;
      vtkWarningMacro(<<"StartTimeStep out of range; setting to " << StartTimeStep);
    }
  }
  else {
    vtkErrorMacro(<<"Input information has no TIME_STEPS set");
    return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkVofAdvect::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
				    vtkInformationVector **inputVector,
				    vtkInformationVector *outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);

  if(this->FirstIteration) {
    TerminationTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    if(this->TerminationTime > this->InputTimeValues.back()) {
      this->TerminationTime = this->InputTimeValues.back();
    }
    TerminationTimeStep = findClosestTimeStep(TerminationTime, InputTimeValues);

    if (TerminationTimeStep < 0) {
      TerminationTimeStep = 0;
      vtkWarningMacro(<<"TerminationTimeStep out of range; setting to " << TerminationTimeStep);
    }
    if (TerminationTimeStep > InputTimeValues.size()-1) {
      TerminationTimeStep = InputTimeValues.size()-1;
      vtkWarningMacro(<<"TerminationTimeStep out of range; setting to " << TerminationTimeStep);
    }
    CurrentTimeStep = StartTimeStep;
  }

  if (StartTimeStep <= TerminationTimeStep) {
    int numInputs = this->GetNumberOfInputPorts();
    for (int i = 0; i < numInputs; i++) {
      vtkInformation *inInfo = inputVector[i]->GetInformationObject(0);

      if (this->CurrentTimeStep < static_cast<int>(this->InputTimeValues.size())) {
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), 
		    this->InputTimeValues[this->CurrentTimeStep]);
      }
    }
  }

  return 1;
}
//----------------------------------------------------------------------------
int vtkVofAdvect::RequestData(vtkInformation *request,
			      vtkInformationVector **inputVector,
			      vtkInformationVector *outputVector)
{
  if (StartTimeStep > TerminationTimeStep) 
    return 1;

  // get the info objects
  vtkInformation *inInfoVelocity = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfoVof = inputVector[1]->GetInformationObject(0);
  vtkInformation *inInfoParticles = inputVector[2]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input
  vtkRectilinearGrid *inputVelocityGrid = vtkRectilinearGrid::
    SafeDownCast(inInfoVelocity->Get(vtkDataObject::DATA_OBJECT()));
  vtkRectilinearGrid *inputVofGrid = vtkRectilinearGrid::
    SafeDownCast(inInfoVof->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *inputParticles = vtkPolyData::
    SafeDownCast(inInfoParticles->Get(vtkDataObject::DATA_OBJECT()));

  int processId = 0;
  int numProcesses = Controller->GetNumberOfProcesses();
  if (numProcesses > 0) { 
    processId = Controller->GetLocalProcessId();
  }
  
  if (this->FirstIteration) {

    vtkDataArray *coords[3];
    getGridCenterCoords(inputVelocityGrid, vtkDataSetAttributes::VECTORS, coords);
    XCoords = coords[0];
    YCoords = coords[1];
    ZCoords = coords[2];

    vtkPoints *seeds = inputParticles->GetPoints();
    Particles.clear();
    for (int i = 0; i < seeds->GetNumberOfPoints(); ++i) {
     
      double p[3];
      seeds->GetPoint(i, p);
      Particles.push_back(make_f3u1_t(p[0], p[1], p[2], i));
    }

    PProcessId = std::vector<short>(Particles.size(), processId);

    //----------------------------------------------------
    // find out which processes are adjacent to this one
    if (numProcesses > 0) {

      inputVelocityGrid->GetBounds(&MyBounds[0]);

      // ------------------------
      std::vector<vtkIdType> RecvLengths(numProcesses);
      std::vector<vtkIdType> RecvOffsets(numProcesses);
      for (int i = 0; i < numProcesses; ++i) {
	RecvLengths[i] = 6;
	RecvOffsets[i] = i*6;
      }

      int MyExtent[6];
      inputVelocityGrid->GetExtent(MyExtent);
      std::vector<int> AllExtents(6*numProcesses);
      Controller->AllGatherV(&MyExtent[0], &AllExtents[0], 6, &RecvLengths[0], &RecvOffsets[0]);

      std::vector<double> AllBounds(6*numProcesses);
      Controller->AllGatherV(&MyBounds[0], &AllBounds[0], 6, &RecvLengths[0], &RecvOffsets[0]);
    
      findGlobalExtents(AllExtents, GlobalExtents);
      findGlobalBounds(AllBounds, GlobalBounds);

      NeighborProcesses.clear();
      NeighborProcesses.resize(6);
      findNeighbors(MyExtent, GlobalExtents, AllExtents, NeighborProcesses);
      // ------------------------

      NumNeighbors = 0;
      for (int i = 0; i < NeighborProcesses.size(); ++i) {
	NumNeighbors += NeighborProcesses[i].size();
      }
    }
    //---------------------------------------------------
  }

  bool finished = this->CurrentTimeStep >= this->TerminationTimeStep;
  if (!finished) {
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
    this->FirstIteration = false;
  }
  else {
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    this->FirstIteration = true;
  }

  //--------------------------------------------------------------------------
  // processing here
  if (StartTimeStep < TerminationTimeStep && 
      CurrentTimeStep < TerminationTimeStep) {

    std::cout << "{" << processId << "} [" << CurrentTimeStep << "] advection start" << std::endl;
    Advect(inputVofGrid, inputVelocityGrid);
    std::cout << "{" << processId << "} [" << CurrentTimeStep << "] advection end" << std::endl;
  }

  //
  if (numProcesses > 0) {

    processId = Controller->GetLocalProcessId();
    // one vector for each side of the process
    std::vector<std::vector<f3u1_t> > particlesToSend(numProcesses);
    // particle process ids hold initial process id of each particle;
    // this is necessary to send labels from advected particles to 
    // the seeds; the local particle id is stored as .id component of
    // the f3u1_t type
    std::vector<std::vector<short> > pidsToSend(numProcesses);
    
    for (int i = 0; i < numProcesses; ++i) {
      particlesToSend[i].resize(0);
      pidsToSend[i].resize(0);
    }

    //------------------------------------------------------------------------
    std::vector<f3u1_t>::iterator it;
    std::vector<f3u1_t> particlesToKeep;
    std::vector<short> pidsToKeep;

    int pidx = 0;
    for (it = Particles.begin(); it != Particles.end(); ++it) {

      int bound = outOfBounds(*it, MyBounds, GlobalBounds);
      if (bound > -1) {
	for (int j = 0; j < NeighborProcesses[bound].size(); ++j) {

	  int neighborId = NeighborProcesses[bound][j];

	  particlesToSend[neighborId].push_back(*it);
	  pidsToSend[neighborId].push_back(PProcessId[pidx]);
	}
      }
      else {
	particlesToKeep.push_back(*it);
	pidsToKeep.push_back(PProcessId[pidx]);
      }

      ++pidx;
    }
    Particles = particlesToKeep;
    PProcessId = pidsToKeep;
    //------------------------------------------------------------------------

    std::vector<f3u1_t> particlesToRecv;
    sendData(particlesToSend, particlesToRecv, numProcesses, Controller);
    std::vector<short> pidsToRecv;
    sendData(pidsToSend, pidsToRecv, numProcesses, Controller);

    // insert the paricles that are within the domain
    for (int i = 0; i < particlesToRecv.size(); ++i) {
      int within = withinBounds(particlesToRecv[i], MyBounds);
      if (within) {
	Particles.push_back(particlesToRecv[i]);
	PProcessId.push_back(pidsToRecv[i]); // _pid
      }
    }
  }

  if(StartTimeStep <= TerminationTimeStep &&
     CurrentTimeStep <= TerminationTimeStep) {
    CurrentTimeStep++;
  }

  //---------------
  // get the output
  if (CurrentTimeStep == TerminationTimeStep) {

    vtkSmartPointer<vtkPoints> output_endPoints = vtkSmartPointer<vtkPoints>::New();
    vtkIntArray *particleId = vtkIntArray::New();
    particleId->SetName("ParticleId");
    particleId->SetNumberOfComponents(1);
    particleId->SetNumberOfTuples(Particles.size()); 
    vtkShortArray *particleProcessId = vtkShortArray::New();
    particleProcessId->SetName("ParticleProcessId");
    particleProcessId->SetNumberOfComponents(1);
    particleProcessId->SetNumberOfTuples(Particles.size());

    for (int i = 0; i < Particles.size(); ++i) {

      double p[3];
      p[0] = Particles[i].x;
      p[1] = Particles[i].y;
      p[2] = Particles[i].z;
      output_endPoints->InsertNextPoint(p);

      particleId->SetValue(i, Particles[i].id);
      particleProcessId->SetValue(i, PProcessId[i]);
    }

    vtkPolyData *output = 
      vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));    
    output->SetPoints(output_endPoints);
    output->GetPointData()->AddArray(particleId);
    output->GetPointData()->AddArray(particleProcessId);
  }  
  return 1;
}

int vtkVofAdvect::Advect(vtkRectilinearGrid *inputVofGrid,
			 vtkRectilinearGrid *inputVelocityGrid)
{
  //------------------------------------------------
  // discard particles that move to cells with f = 0
  vtkDataArray *vofArray;
  getGridData(inputVofGrid, vtkDataSetAttributes::SCALARS, &vofArray);
  int res[3] = {XCoords->GetNumberOfTuples(),
		YCoords->GetNumberOfTuples(),
		ZCoords->GetNumberOfTuples()};
  
  if (vofArray->GetDataType() == VTK_FLOAT) {
    discardStrayParticles(vtkFloatArray::SafeDownCast(vofArray)->GetPointer(0), res, 
    			  vtkFloatArray::SafeDownCast(XCoords)->GetPointer(0), 
    			  vtkFloatArray::SafeDownCast(YCoords)->GetPointer(0), 
    			  vtkFloatArray::SafeDownCast(ZCoords)->GetPointer(0), Particles);
  }
  else if (vofArray->GetDataType() == VTK_DOUBLE) {
    discardStrayParticles(vtkDoubleArray::SafeDownCast(vofArray)->GetPointer(0), res, 
    			  vtkDoubleArray::SafeDownCast(XCoords)->GetPointer(0), 
    			  vtkDoubleArray::SafeDownCast(YCoords)->GetPointer(0), 
    			  vtkDoubleArray::SafeDownCast(ZCoords)->GetPointer(0), Particles);
  }

  vtkDataArray *velocityArray;
  getGridData(inputVelocityGrid, vtkDataSetAttributes::VECTORS, &velocityArray);  

  float dt = InputTimeValues[CurrentTimeStep+1] - InputTimeValues[CurrentTimeStep];
  // Time step in the simulation data might be incorrent, in which case we set it
  // manually
  if (TimeStepDelta != 0.0) {
    dt = TimeStepDelta;    
  }

  if (velocityArray->GetDataType() == VTK_FLOAT) {    
    advectParticles(vtkFloatArray::SafeDownCast(velocityArray)->GetPointer(0), res, 
  		    vtkFloatArray::SafeDownCast(XCoords)->GetPointer(0), 
  		    vtkFloatArray::SafeDownCast(YCoords)->GetPointer(0), 
  		    vtkFloatArray::SafeDownCast(ZCoords)->GetPointer(0), 
  		    dt, Particles);
  }
  else if (velocityArray->GetDataType() == VTK_DOUBLE) {
    advectParticles(vtkDoubleArray::SafeDownCast(velocityArray)->GetPointer(0), res, 
  		    vtkDoubleArray::SafeDownCast(XCoords)->GetPointer(0), 
  		    vtkDoubleArray::SafeDownCast(YCoords)->GetPointer(0), 
  		    vtkDoubleArray::SafeDownCast(ZCoords)->GetPointer(0), 
  		    dt, Particles);
  }

  return 1;
}

////////// External Operators /////////////

void vtkVofAdvect::PrintSelf(ostream &os, vtkIndent indent)
{
}
