#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkMPIController.h"

#include "vtkVofComponents.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <map>

#include "helper_math.h"

vtkStandardNewMacro(vtkVofComponents);

namespace
{
  int uf_root(std::vector<int> &id, int i)
  {
    while (i != id[i]) {
      id[i] = id[id[i]];
      i = id[i];
    }
    return i;
  }

  bool uf_find(std::vector<int> &id, int p, int q)
  {
    return uf_root(id, p) == uf_root(id, q);
  }

  void uf_unite(std::vector<int> &id, int p, int q)
  {
    int i = uf_root(id, p);
    int j = uf_root(id, q);
    id[i] = j;
  }

  // multiprocess
void findGlobalExtent(std::vector<int> &allExtents, 
		      int globalExtent[6])
{
  globalExtent[0] = globalExtent[2] = globalExtent[4] = std::numeric_limits<int>::max();
  globalExtent[1] = globalExtent[3] = globalExtent[5] = - globalExtent[0];

  for (int i = 0; i < allExtents.size()/6; ++i) {
    if (globalExtent[0] > allExtents[i*6+0]) globalExtent[0] = allExtents[i*6+0];
    if (globalExtent[1] < allExtents[i*6+1]) globalExtent[1] = allExtents[i*6+1];
    if (globalExtent[2] > allExtents[i*6+2]) globalExtent[2] = allExtents[i*6+2];
    if (globalExtent[3] < allExtents[i*6+3]) globalExtent[3] = allExtents[i*6+3];
    if (globalExtent[4] > allExtents[i*6+4]) globalExtent[4] = allExtents[i*6+4];
    if (globalExtent[5] < allExtents[i*6+5]) globalExtent[5] = allExtents[i*6+5];
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

void findNeighbors(const int localExtents[6], 
		   const int globalExtents[6], 
		   const std::vector<int> &allExtents,
		   std::vector<std::vector<int> > &neighbors)
{
  const int numDims = 3;
  const int numSides = 6;
  
  for (int i = 0; i < numDims; ++i) {

    if (localExtents[i*2+0] > globalExtents[i*2+0]) {
      for (int j = 0; j < allExtents.size()/numSides; ++j) {
	
	if (localExtents[i*2+0] <= allExtents[j*numSides+i*2+1] &&
	    localExtents[i*2+1] > allExtents[j*numSides+i*2+1] &&
	    localExtents[((i+1)%3)*2+0] < allExtents[j*numSides+((i+1)%3)*2+1] &&
	    localExtents[((i+1)%3)*2+1] > allExtents[j*numSides+((i+1)%3)*2+0] &&
	    localExtents[((i+2)%3)*2+0] < allExtents[j*numSides+((i+2)%3)*2+1] &&
	    localExtents[((i+2)%3)*2+1] > allExtents[j*numSides+((i+2)%3)*2+0]) {

	  neighbors[i*2+0].push_back(j);
	}
      }
    }
    if (localExtents[i*2+1] < globalExtents[i*2+1]) { 
      for (int j = 0; j < allExtents.size()/numSides; ++j) {

	if (localExtents[i*2+1] >= allExtents[j*numSides+i*2+0] &&
	    localExtents[i*2+0] < allExtents[j*numSides+i*2+0] &&
	    localExtents[((i+1)%3)*2+0] < allExtents[j*numSides+((i+1)%3)*2+1] &&
	    localExtents[((i+1)%3)*2+1] > allExtents[j*numSides+((i+1)%3)*2+0] &&
	    localExtents[((i+2)%3)*2+0] < allExtents[j*numSides+((i+2)%3)*2+1] &&
	    localExtents[((i+2)%3)*2+1] > allExtents[j*numSides+((i+2)%3)*2+0]) {

	  neighbors[i*2+1].push_back(j);
	}
      }
    }
  }
}

}

static const double g_emf0 = 0.000001;
static int g_res[3];
static int g_labelId = 0;

template<typename T>
static void grow(std::vector<T> &data, const int &idx,
		 std::vector<float> &labels)
{
  T tmp = data[idx];

  if (tmp <= g_emf0) {
    return;
  }

  labels[idx] = g_labelId;
  data[idx] = 0.0;

  int i = idx%g_res[0];
  int j = (idx/g_res[0])%g_res[1];
  int k = idx/(g_res[0]*g_res[1]);

  if (i + 1 < g_res[0]) {
    grow(data, idx+1, labels);
  }
  if (i > 0) {
    grow(data, idx-1, labels);
  }
  if (j + 1 < g_res[1]) {
    grow(data, idx+g_res[0], labels);
  }
  if (j > 0) {
    grow(data, idx-g_res[0], labels);
  }
  if (k + 1 < g_res[2]) {
    grow(data, idx+g_res[0]*g_res[1], labels);
  }
  if (k > 0) {
    grow(data, idx-g_res[0]*g_res[1], labels);
  }
}

template<typename T>
static bool extractComponents_pred(T f)
{
  return f > g_emf0;
}

template<typename T>
static void extractComponents(const T *vofField,
			      const int res[3],
			      float *labelField)
{
  g_res[0] = res[0];
  g_res[1] = res[1];
  g_res[2] = res[2];

  int numCells = res[0]*res[1]*res[2];
  std::vector<T> fieldTmp(vofField, vofField+numCells);
  std::vector<float> labelFieldTmp(numCells, -1.0f);

  int i = 0;
  while (i < numCells) {

    typename std::vector<T>::iterator it =
      std::find_if(fieldTmp.begin()+i, fieldTmp.end(),
		   extractComponents_pred<T>);

    if (it == fieldTmp.end()) {
      break;
    }

    int idx = it - fieldTmp.begin();

    if (idx > i) {
      i = idx;
    }
    else {
      ++i;
    }
    grow(fieldTmp, idx, labelFieldTmp);
    ++g_labelId;
  }
  std::copy(labelFieldTmp.begin(), labelFieldTmp.end(), labelField);
}

static
void unifyLabelsInProcess(std::vector<std::vector<int> > &NeighborProcesses,
			  const int myExtent[6], int cellRes[3], vtkFloatArray *labels,
			  std::vector<std::vector<float4> > &labelsToRecv,
			  std::vector<int> &labelOffsets, int processId,
			  std::vector<int> &allLabels)
{
  const int NUM_SIDES = 6;
  int nidx = 0;
  for (int i = 0; i < NUM_SIDES; ++i) {
    for (int j = 0; j < NeighborProcesses[i].size(); ++j) {
      for (int s = 0; s < labelsToRecv[nidx].size(); ++s) {

	int x = labelsToRecv[nidx][s].x;
	int y = labelsToRecv[nidx][s].y;
	int z = labelsToRecv[nidx][s].z;

	if (x >= myExtent[0] && x <= myExtent[1] &&
	    y >= myExtent[2] && y <= myExtent[3] &&
	    z >= myExtent[4] && z <= myExtent[5]) {
	    
	  x -= myExtent[0];
	  y -= myExtent[2];
	  z -= myExtent[4];

	  int idx = x + y*cellRes[0] + z*cellRes[0]*cellRes[1];
	  int myLabel = labels->GetValue(idx) + labelOffsets[processId];

	  int neighborLabel = labelsToRecv[nidx][s].w + labelOffsets[NeighborProcesses[i][j]];

	  if (labels->GetValue(idx) > -1) {
	    if (!uf_find(allLabels, myLabel, neighborLabel)) {
	      uf_unite(allLabels, myLabel, neighborLabel);
	    }
	  }
	}
      }
      ++nidx;
    }
  }
}

static
void unifyLabelsInDomain(std::vector<int> &allLabelUnions, int numAllLabels,
			 std::vector<int> &allLabels, vtkFloatArray *labels,
			 std::vector<int> &labelOffsets, int processId)
{
  for (int i = 0; i < allLabelUnions.size(); ++i) {

    int labelId = i%numAllLabels;
    if (allLabelUnions[i] != labelId) {
      if (!uf_find(allLabels, allLabelUnions[i], labelId)) {
	uf_unite(allLabels, allLabelUnions[i], labelId);
      }
    }
  }

  for (int i = 0; i < allLabels.size(); ++i) {
    if (allLabels[i] != i) {
      int rootId = uf_root(allLabels, i);
      allLabels[i] = rootId;
    }
  }
  std::map<int,int> labelMap;
  int labelId = 0;
  for (int i = 0; i < allLabels.size(); ++i) {
    if (labelMap.find(allLabels[i]) == labelMap.end()) {
      labelMap[allLabels[i]] = labelId;
      ++labelId;
    }
  }

  float *labels_ptr = labels->GetPointer(0);
  for (int i = 0; i < labels->GetNumberOfTuples(); ++i) {
    if (labels_ptr[i] > -1) {
      int label = labels_ptr[i] + labelOffsets[processId];
      label = allLabels[label];
      labels_ptr[i] = labelMap[label];
    }
  }
}

static void prepareLabelsToSend(std::vector<std::vector<int> > &NeighborProcesses,
				const int myExtent[6], int cellRes[3], vtkFloatArray *labels,
				std::vector<std::vector<float4> > &labelsToSend)
{
  const int NUM_SIDES = 6;
  const int slabs[NUM_SIDES][6] =
    {{0,1,                       0,cellRes[1]-1,            0,cellRes[2]-1},
     {cellRes[0]-2,cellRes[0]-1, 0,cellRes[1]-1,            0,cellRes[2]-1},
     {0,cellRes[0]-1,            0,1,                       0,cellRes[2]-1},
     {0,cellRes[0]-1,            cellRes[1]-2,cellRes[1]-1, 0,cellRes[2]-1},
     {0,cellRes[0]-1,            0,cellRes[1]-1,            0,1},
     {0,cellRes[0]-1,            0,cellRes[1]-1,            cellRes[2]-2,cellRes[2]-1}};

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
}

//-----------------------------------------------------------------------------
vtkVofComponents::vtkVofComponents() :
  NumGhostLevels(1)
{
  Controller = vtkMPIController::New();
}

//-----------------------------------------------------------------------------
vtkVofComponents::~vtkVofComponents()
{
  Controller->Delete();
}

//----------------------------------------------------------------------------
void vtkVofComponents::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}


//----------------------------------------------------------------------------
void vtkVofComponents::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
int vtkVofComponents::FillInputPortInformation( int port, vtkInformation* info )
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
      return 0;
    }
  if ( port == 0 )
    {
      info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
      return 1;
    }
  return 0;
}

void vtkVofComponents::GetGlobalContext(vtkInformation *inInfo)
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
}

void vtkVofComponents::ExtractComponents(vtkRectilinearGrid *vof,
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
int vtkVofComponents::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
					  vtkInformationVector **inputVector,
					  vtkInformationVector *outputVector)
{
  // set one ghost level -----------------------------------------------------
  int numInputs = this->GetNumberOfInputPorts();
  for (int i = 0; i < numInputs; i++) {
    vtkInformation *inInfo = inputVector[i]->GetInformationObject(0);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), NumGhostLevels);
  }
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), NumGhostLevels);
}
//----------------------------------------------------------------------------
int vtkVofComponents::RequestData(vtkInformation *vtkNotUsed(request),
				  vtkInformationVector **inputVector,
				  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  if (Controller->GetCommunicator() != 0) {
    GetGlobalContext(inInfo);
  }

  // get the input and output
  vtkRectilinearGrid *input = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkRectilinearGrid *output = vtkRectilinearGrid::
    SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    
  ExtractComponents(input, output);

  return 1;
}

////////// External Operators /////////////

void vtkVofComponents::PrintSelf(ostream &os, vtkIndent indent)
{
}
