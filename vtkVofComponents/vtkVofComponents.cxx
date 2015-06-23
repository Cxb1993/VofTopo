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

  void findNeighbors(const double myBounds[6],
		     const double globalBounds[6],
		     const std::vector<double> &allBounds,
		     std::vector<std::vector<int> > &neighbors)
  {
    const int numDims = 3;
    const int numSides = 6;

    for (int i = 0; i < numDims; ++i) {

      if (myBounds[i*2+0] > globalBounds[i*2+0]) {
	for (int j = 0; j < allBounds.size()/numSides; ++j) {

	  if (myBounds[i*2+0] == allBounds[j*numSides+i*2+1] &&
	      myBounds[((i+1)%3)*2+0] < allBounds[j*numSides+((i+1)%3)*2+1] &&
	      myBounds[((i+1)%3)*2+1] > allBounds[j*numSides+((i+1)%3)*2+0] &&
	      myBounds[((i+2)%3)*2+0] < allBounds[j*numSides+((i+2)%3)*2+1] &&
	      myBounds[((i+2)%3)*2+1] > allBounds[j*numSides+((i+2)%3)*2+0]) {

	    neighbors[i*2+0].push_back(j);
	  }
	}
      }
      if (myBounds[i*2+1] < globalBounds[i*2+1]) {
	for (int j = 0; j < allBounds.size()/numSides; ++j) {

	  if (myBounds[i*2+1] == allBounds[j*numSides+i*2+0] &&
	      myBounds[((i+1)%3)*2+0] < allBounds[j*numSides+((i+1)%3)*2+1] &&
	      myBounds[((i+1)%3)*2+1] > allBounds[j*numSides+((i+1)%3)*2+0] &&
	      myBounds[((i+2)%3)*2+0] < allBounds[j*numSides+((i+2)%3)*2+1] &&
	      myBounds[((i+2)%3)*2+1] > allBounds[j*numSides+((i+2)%3)*2+0]) {

	    neighbors[i*2+1].push_back(j);
	  }
	}
      }
    }
  }
  
  bool adjointPoint(int x, int y, int z, int extent[6]) 
  {
    if (x >= extent[0] && x <= extent[1] && 
	y >= extent[2] && y <= extent[3] && 
	z >= extent[4] && z <= extent[5]) {
      return true;
    }
    return false;
  }
  bool adjointCell(int x, int y, int z, int extent[6], int side,
		   int &mx, int &my, int &mz) 
  {
    if (side == 0) { // left neighbor
      if (x == extent[0]-1 && 
	  y >= extent[2] && y < extent[3] && 
	  z >= extent[4] && z < extent[5]) {
	mx = 0;
	my = y - extent[2];
	mz = z - extent[4];
	return true;
      }
    }
    else if (side == 1) { // right corner
      if (x == extent[1] && 
	  y >= extent[2] && y < extent[3] && 
	  z >= extent[4] && z < extent[5]) {
	mx = extent[1]-1;
	my = y - extent[2];
	mz = z - extent[4];
	return true;
      }
    }
    if (side == 2) { // bottom neighbor
      if (x >= extent[0] && x < extent[1] &&
	  y == extent[2]-1 && 
	  z >= extent[4] && z < extent[5]) {
	mx = x - extent[0];
	my = 0;
	mz = z - extent[4];
	return true;
      }
    }
    else if (side == 3) { // top corner
      if (x >= extent[0] && x < extent[1] &&
	  y == extent[3] && 
	  z >= extent[4] && z < extent[5]) {
	mx = x - extent[0];
	my = extent[3]-1;
	mz = z - extent[4];
	return true;
      }
    }
    if (side == 4) { // back neighbor
      if (x >= extent[0] && x < extent[1] &&
	  y >= extent[2] && y < extent[3] && 
	  z == extent[4]-1) {
	mx = x - extent[0];
	my = y - extent[2];
	mz = 0;
	return true;
      }
    }
    else if (side == 5) { // front corner
      if (x >= extent[0] && x < extent[1] &&
	  y >= extent[2] && y < extent[3] && 
	  z == extent[5]) {
	mx = x - extent[0];
	my = y - extent[2];
	mz = extent[5]-1;
	return true;
      }
    }
    return false;
  }
}

//-----------------------------------------------------------------------------
vtkVofComponents::vtkVofComponents()
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

int vtkVofComponents::RequestInformation(vtkInformation *vtkNotUsed(request),
					 vtkInformationVector **inputVector,
					 vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

  int extent[6];
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

  return 1;
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

//----------------------------------------------------------------------------
int vtkVofComponents::RequestData(vtkInformation *vtkNotUsed(request),
				  vtkInformationVector **inputVector,
				  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkRectilinearGrid *input = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  int inputRes[3];
  input->GetDimensions(inputRes);
  int outputRes[3] = {inputRes[0], inputRes[1], inputRes[2]};

  // determine if data is cell-based or point-based
  vtkDataArray *pointData = input->GetPointData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  vtkDataArray *cellData = input->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  vtkDataArray *data;

  // cell data
  if (pointData == 0 && cellData != 0) {
    data = cellData;
    outputRes[0] -= 1;
    outputRes[1] -= 1;
    outputRes[2] -= 1;
  }
  // point data
  else if (pointData != 0 && cellData == 0) {
    data = pointData;
  }

  vtkFloatArray *labels = vtkFloatArray::New();
  labels->SetName("Labels");
  labels->SetNumberOfComponents(1);
  labels->SetNumberOfTuples(outputRes[0]*outputRes[1]*outputRes[2]);
  for (int i = 0; i < labels->GetNumberOfTuples(); ++i) {
    labels->SetValue(i, -1.0f);
  }

  g_labelId = 0;
  // determine if data is float or double
  if (data->IsA("vtkFloatArray")) {
    extractComponents(vtkFloatArray::SafeDownCast(data)->GetPointer(0),
  		      outputRes, labels->GetPointer(0));
  }
  else if (data->IsA("vtkDoubleArray")) {
    extractComponents(vtkDoubleArray::SafeDownCast(data)->GetPointer(0),
  		      outputRes, labels->GetPointer(0));
  }

  //--------------------------------------------------------------------------
  // send number of labels to other processes
  int numProcesses = this->Controller->GetNumberOfProcesses();

  if (numProcesses > 0) { // if running in mpi

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

    // -----------------------------------------------------------------------
    // find neighboring processes
    double myBounds[6];
    input->GetBounds(&myBounds[0]);
    std::vector<double> allBounds(numProcesses*6);
    for (int i = 0; i < numProcesses; ++i) {
      recvLengths[i] = 6;
      recvOffsets[i] = i*6;
    }
    Controller->AllGatherV(&myBounds[0], &allBounds[0], 6, &recvLengths[0], &recvOffsets[0]);

    double globalBounds[6];
    findGlobalBounds(allBounds, globalBounds);

    std::vector<std::vector<int> > neighborProcesses;
    neighborProcesses.clear();
    neighborProcesses.resize(6);
    findNeighbors(myBounds, globalBounds, allBounds, neighborProcesses);

    int numNeighbors = 0;
    for (int i = 0; i < neighborProcesses.size(); ++i) {
      numNeighbors += neighborProcesses[i].size();
    }

    // -----------------------------------------------------------------------
    // prepare labelled points to send to neighbors
    int myExtent[6];
    input->GetExtent(myExtent);

    int res[3];
    input->GetDimensions(res);
    if (pointData == 0 && cellData != 0) {
      res[0] -= 1;
      res[1] -= 1;
      res[2] -= 1;
    }

    int lext[6][6] = {{0,0,                0,res[1]-1,         0,res[2]-1},
		      {res[0]-1,res[0]-1,  0,res[1]-1,         0,res[2]-1},
		      {0,res[0]-1,         0,0,                0,res[2]-1},
		      {0,res[0]-1,         res[1]-1,res[1]-1,  0,res[2]-1},
		      {0,res[0]-1,         0,res[1]-1,         0,0},
		      {0,res[0]-1,         0,res[1]-1,         res[2]-1,res[2]-1}};

    std::vector<std::vector<float4> > labelsToSend(6);
    for (int p = 0; p < neighborProcesses.size(); ++p) {
      if (neighborProcesses[p].size() > 0) {

	for (int k = lext[p][4]; k <= lext[p][5]; ++k) {
	  for (int j = lext[p][2]; j <= lext[p][3]; ++j) {
	    for (int i = lext[p][0]; i <= lext[p][1]; ++i) {

	      int idx = i + j*res[0] + k*res[0]*res[1];
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

      for (int j = 0; j < neighborProcesses[i].size(); ++j) {

	const int SEND_LABELS_TAG = 100+processId;
	Controller->Send(&numLabelsToSend[i], 1,
			 neighborProcesses[i][j], SEND_LABELS_TAG);
      }
    }

    // -----------------------------------------------------------------------
    // receive header
    std::vector<int> numLabelsToRecv(numNeighbors);
    int nidx = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < neighborProcesses[i].size(); ++j) {
	numLabelsToRecv[nidx] = 0;
	const int RECV_LABELS_TAG = 100+neighborProcesses[i][j];
	Controller->Receive(&numLabelsToRecv[nidx], 1,
			    neighborProcesses[i][j], RECV_LABELS_TAG);
	++nidx;
      }
    }

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------

    // send the labels to each side
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < neighborProcesses[i].size(); ++j) {
    	const int SEND_LABEL_DATA_TAG = 100+processId;

    	vtkMPICommunicator::Request req;
    	Controller->NoBlockSend((char*)&(labelsToSend[i][0]),
    				labelsToSend[i].size()*sizeof(float4),
    				neighborProcesses[i][j], SEND_LABEL_DATA_TAG, req);
      }
    }
    // -----------------------------------------------------------------------
    // allocate buffers to receive labels from each neighbor
    std::vector<std::vector<float4> > labelsToRecv(numNeighbors);
    for (int i = 0; i < numNeighbors; ++i) {
      labelsToRecv[i].resize(numLabelsToRecv[i]);
    }
    // -----------------------------------------------------------------------
    // receive labels from each neighbor
    vtkMPICommunicator::Request *reqs = new vtkMPICommunicator::Request[numNeighbors];
    nidx = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < neighborProcesses[i].size(); ++j) {
    	const int RECV_LABEL_DATA_TAG = 100+neighborProcesses[i][j];

    	Controller->NoBlockReceive((char*)&(labelsToRecv[nidx][0]),
    				   labelsToRecv[nidx].size()*sizeof(float4),
    				   neighborProcesses[i][j], RECV_LABEL_DATA_TAG, reqs[nidx]);
    	++nidx;
      }
    }
    Controller->WaitAll(numNeighbors, reqs);

    std::vector<int> allLabels(numAllLabels);
    for (int i = 0; i < allLabels.size(); ++i) {
      allLabels[i] = i;
    }

    nidx = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < neighborProcesses[i].size(); ++j) {
    	for (int s = 0; s < labelsToRecv[nidx].size(); ++s) {

    	  int x = labelsToRecv[nidx][s].x;
    	  int y = labelsToRecv[nidx][s].y;
    	  int z = labelsToRecv[nidx][s].z;
    	  int neighborLabel = labelsToRecv[nidx][s].w + labelOffsets[neighborProcesses[i][j]];
	  
	  int mx, my, mz;

	  if (pointData != 0 && cellData == 0) {
	    if (!adjointPoint(x, y, z, myExtent)) {
	      continue;
	    }
	    mx = x - myExtent[0];
	    my = y - myExtent[2];
	    mz = z - myExtent[4];
	  }
	  if (pointData == 0 && cellData != 0) {
	    if (!adjointCell(x, y, z, myExtent, i, mx, my, mz)) {
	      continue;
	    }	    
	  }

    	  int idx = mx + my*res[0] + mz*res[0]*res[1];
    	  int myLabel = labels->GetValue(idx) + labelOffsets[processId];

	  if (labels->GetValue(idx) > -1) {
	    if (!uf_find(allLabels, myLabel, neighborLabel)) {
	      uf_unite(allLabels, myLabel, neighborLabel);
	    }
	  }
    	}
    	++nidx;
      }
    }
    for (int i = 0; i < numProcesses; ++i) {
      recvLengths[i] = numAllLabels;
      recvOffsets[i] = i*numAllLabels;
    }
    std::vector<int> allLabelUnions(numAllLabels*numProcesses);
    Controller->AllGatherV(&allLabels[0], &allLabelUnions[0], numAllLabels, 
    			   &recvLengths[0], &recvOffsets[0]);

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

  vtkRectilinearGrid *output = vtkRectilinearGrid::
    SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  output->SetExtent(input->GetExtent());
  output->SetXCoordinates(input->GetXCoordinates());
  output->SetYCoordinates(input->GetYCoordinates());
  output->SetZCoordinates(input->GetZCoordinates());

  if (pointData == 0 && cellData != 0) {
    output->GetCellData()->AddArray(labels);
    output->GetCellData()->SetActiveScalars("Labels");
  }
  if (pointData != 0 && cellData == 0) {
    output->GetPointData()->AddArray(labels);
    output->GetPointData()->SetActiveScalars("Labels");
  }

  return 1;
}

////////// External Operators /////////////

void vtkVofComponents::PrintSelf(ostream &os, vtkIndent indent)
{
}
