#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkPolyData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkShortArray.h"
#include "vtkCharArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPoints.h"
#include "vtkMPIController.h"

#include "vtkVofLabelPoints.h"

#include <vector>
#include <cmath>

vtkStandardNewMacro(vtkVofLabelPoints);

namespace 
{
// unsafe! does not take cell sizes into account
  int findClosestCoordinate(vtkDataArray *coords, double p) 
  {
    const int numCoordinates = coords->GetNumberOfTuples(); 
    int coord = 0;
    double dist = std::abs(coords->GetComponent(0,0)-p);
    for (int i = 1; i < numCoordinates; ++i) {

      double curr_dist = std::abs(coords->GetComponent(i,0)-p);
      if (curr_dist < dist) { 
	dist = curr_dist;
	coord = i;
      }
    }
    return coord;
  }

  int findCell(vtkDataArray *coords, double p) 
  {
    const int numCoordinates = coords->GetNumberOfTuples(); 
    int coord = 0;
    for (int i = 0; i < numCoordinates-1; ++i) {
      if (p >= coords->GetComponent(i,0) && p <= coords->GetComponent(i+1,0)) {
	coord = i;
	break;
      }
    }
    return coord;
  }

  typedef struct {
    int particleId;
    float label;
  } i1f1_t;
}

//-----------------------------------------------------------------------------
vtkVofLabelPoints::vtkVofLabelPoints()
{
  this->SetNumberOfInputPorts(3);
  this->Controller = vtkMPIController::New();
}

//-----------------------------------------------------------------------------
vtkVofLabelPoints::~vtkVofLabelPoints()
{
}

//----------------------------------------------------------------------------
void vtkVofLabelPoints::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}

//----------------------------------------------------------------------------
void vtkVofLabelPoints::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
int vtkVofLabelPoints::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
  }
  if (port == 1) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  if (port == 2) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
int vtkVofLabelPoints::RequestData(vtkInformation *request,
				   vtkInformationVector **inputVector,
				   vtkInformationVector *outputVector)
{
  int processId = 0;
  int numProcesses = Controller->GetNumberOfProcesses();
  if (numProcesses > 0) { 
    processId = Controller->GetLocalProcessId();
  }

  // get the info objects
  vtkInformation *inInfoComponents = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfoSeeds = inputVector[1]->GetInformationObject(0);
  vtkInformation *inInfoAdvectedParticles = inputVector[2]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input
  vtkRectilinearGrid *inputComponents = vtkRectilinearGrid::
    SafeDownCast(inInfoComponents->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *inputSeeds = vtkPolyData::
    SafeDownCast(inInfoSeeds->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *inputAdvectedParticles = vtkPolyData::
    SafeDownCast(inInfoAdvectedParticles->Get(vtkDataObject::DATA_OBJECT()));

  if (inputComponents == 0 || inputSeeds == 0 || inputAdvectedParticles == 0) {
    vtkErrorMacro("One of the inputs is empty");
    return 0;
  }

  vtkCharArray *ifacePoints = vtkCharArray::
    SafeDownCast(inputSeeds->GetPointData()->GetArray("InterfacePoints"));
  if (ifacePoints != NULL) std::cout << "BANGLA" << std::endl;
  else std::cout << "NIE BANGLA" << std::endl;

  // output labels for seed points -------------------------------------------
  vtkFloatArray *seedPointLabels;
  float *seedPointLabels_ptr = 0;
  int numSeeds = inputSeeds->GetNumberOfPoints();
  // if (numSeeds > 0) {
    seedPointLabels = vtkFloatArray::New();
    seedPointLabels->SetName("Labels");
    seedPointLabels->SetNumberOfComponents(1);
    seedPointLabels->SetNumberOfTuples(numSeeds);
    seedPointLabels_ptr = seedPointLabels->GetPointer(0);

    for (int i = 0; i < numSeeds; ++i) {
      seedPointLabels_ptr[i] = -1.0f;
    }
  // }

  // output labels for advected particles ------------------------------------
  vtkFloatArray *advectedParticleLabels;
  float *advectedParticleLabels_ptr = 0;
  int numAdvectedParticles = inputAdvectedParticles->GetNumberOfPoints();
  // if (numAdvectedParticles > 0) {
    advectedParticleLabels = vtkFloatArray::New();
    advectedParticleLabels->SetName("Labels");
    advectedParticleLabels->SetNumberOfComponents(1);
    advectedParticleLabels->SetNumberOfTuples(numAdvectedParticles);
    advectedParticleLabels_ptr = advectedParticleLabels->GetPointer(0);
  // }

  std::vector<std::vector<i1f1_t> > labelsToSend;
  labelsToSend.resize(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    labelsToSend[i].resize(0);
  }


  if (numAdvectedParticles > 0) {

    int *particleId_ptr = vtkIntArray::
      SafeDownCast(inputAdvectedParticles->GetPointData()->GetArray("ParticleId"))->GetPointer(0);
    short *particleProcessId_ptr = vtkShortArray::
      SafeDownCast(inputAdvectedParticles->GetPointData()->GetArray("ParticleProcessId"))->
      GetPointer(0);

    if (particleId_ptr == 0 || particleProcessId_ptr == 0) {
      vtkErrorMacro("One fo the Id arrays of advected particles is empty");
      return 0;
    }

    //------------------------------------------------------------------------
    vtkDataArray *pointData = inputComponents->GetPointData()->GetArray("Labels");
    vtkDataArray *cellData = inputComponents->GetCellData()->GetArray("Labels");
    vtkFloatArray *data;

    if (pointData == 0 && cellData != 0) { // cell data
      data = vtkFloatArray::SafeDownCast(cellData);
    }
    else if (pointData != 0 && cellData == 0) { // point data
      data = vtkFloatArray::SafeDownCast(pointData);
    }
    else {
      vtkErrorMacro("Can't determine if data is point-based or cell-based");
      return 0;
    }

    float *componentLabels_ptr = data->GetPointer(0);
    if (componentLabels_ptr == 0) {
      vtkErrorMacro("Component labels array is empty");
      return 0;
    }
    vtkDataArray *coords[3] = {inputComponents->GetXCoordinates(),
			       inputComponents->GetYCoordinates(),
			       inputComponents->GetZCoordinates()};
    int res[3];
    inputComponents->GetDimensions(res);
    if (pointData == 0 && cellData != 0) { // cell data
      res[0] -= 1;
      res[1] -= 1;
      res[2] -= 1;
    }

    vtkPoints *advectedPoints = inputAdvectedParticles->GetPoints();    
    for (int i = 0; i < numAdvectedParticles; ++i) {

      double p[3];
      advectedPoints->GetPoint(i, p);
      int coord[3];
      for (int j = 0; j < 3; ++j) {
	if (pointData != 0 && cellData == 0) {
	  coord[j] = findClosestCoordinate(coords[j], p[j]);
	}
	else if (pointData == 0 && cellData != 0) {
	  coord[j] = findCell(coords[j], p[j]);  
	}
      }
      int idx = coord[0] + coord[1]*res[0] + coord[2]*res[0]*res[1];
      // component label for the given particle
      advectedParticleLabels_ptr[i] = componentLabels_ptr[idx];

      int particleId = particleId_ptr[i];
      if (particleId < 0) { 
	particleId = -1*particleId - 1;
      }

      if (numProcesses > 0) {

	int particleProcId = particleProcessId_ptr[i];
	// if the particle was seeded in this process, assign the label to corresponding seed
	if (particleProcId == processId) {
	  seedPointLabels_ptr[particleId] = componentLabels_ptr[idx];
	}
	else { // if the particle was seeded in other process, store the label
	  i1f1_t labelId = {particleId, componentLabels_ptr[idx]};
	  labelsToSend[particleProcId].push_back(labelId);
	}
      }
      else {
	seedPointLabels_ptr[particleId] = componentLabels_ptr[idx];
      }
    }
  }

  //--------------------------------------------------------------------------
  if (numProcesses > 0) {
    // send labels to particle seeds
    std::vector<int> numLabelsToSend(numProcesses);
    std::vector<int> numLabelsToRecv(numProcesses);
    std::vector<i1f1_t> allLabelsToSend;
    allLabelsToSend.resize(0);
    for (int i = 0; i < numProcesses; ++i) {
      numLabelsToSend[i] = labelsToSend[i].size();
      numLabelsToRecv[i] = 0;
      
      for (int j = 0; j < labelsToSend[i].size(); ++j) {
	allLabelsToSend.push_back(labelsToSend[i][j]);
      }
    }

    std::vector<int> RecvLengths(numProcesses);
    std::vector<int> RecvOffsets(numProcesses);
    int numAllLabelsToRecv = 0;
    for (int i = 0; i < numProcesses; ++i) {
      Controller->Scatter((int*)&numLabelsToSend[0], (int*)&numLabelsToRecv[i], 1, i);

      RecvOffsets[i] = numAllLabelsToRecv;
      RecvLengths[i] = numLabelsToRecv[i]*sizeof(i1f1_t);
      numAllLabelsToRecv += numLabelsToRecv[i];
    }

    std::vector<i1f1_t> labelsToRecv(numAllLabelsToRecv);

    std::vector<vtkIdType> SendLengths(numProcesses);
    std::vector<vtkIdType> SendOffsets(numProcesses);
    int offset = 0;
    for (int i = 0; i < numProcesses; ++i) {
      SendLengths[i] = numLabelsToSend[i]*sizeof(i1f1_t);
      SendOffsets[i] = offset;
      offset += numLabelsToSend[i]*sizeof(i1f1_t);
    }

    for (int i = 0; i < numProcesses; ++i) {
      Controller->ScatterV((char*)&allLabelsToSend[0], (char*)&labelsToRecv[RecvOffsets[i]], 
			   &SendLengths[0], &SendOffsets[0], RecvLengths[i], i);
    }

    for (int i = 0; i < labelsToRecv.size(); ++i) {
      int particleId = labelsToRecv[i].particleId;
      float label = labelsToRecv[i].label;
      seedPointLabels_ptr[particleId] = label;
    }
  }
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  vtkMultiBlockDataSet *output = 
    vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->SetNumberOfBlocks(2);

  // set output blocks
  int output_startParticles_Id = 0;
  {
    if (numSeeds > 0) {
      inputSeeds->GetPointData()->AddArray(seedPointLabels); 
    }
    output->SetBlock(output_startParticles_Id, inputSeeds);
  }
  int output_endParticles_Id = 1;
  {
    if (numAdvectedParticles > 0) {
      inputAdvectedParticles->GetPointData()->AddArray(advectedParticleLabels);
    }
    output->SetBlock(output_endParticles_Id, inputAdvectedParticles);
  }
  
  return 1;
}

////////// External Operators /////////////

void vtkVofLabelPoints::PrintSelf(ostream &os, vtkIndent indent)
{
}
