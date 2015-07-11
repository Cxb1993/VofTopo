// TODO: for vof values of nodes check if ghost nodes are taken 
// into account in grid dimensions
// update: ghost cells lead to many problems, so I'll avoid them for now

#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSetAttributes.h" 
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkShortArray.h"
#include "vtkPointSet.h"
#include "vtkStreamingDemandDrivenPipeline.h"
// #include "vtkMPIController.h"

#include "vtkVofSeedPoints.h"

#include <iostream>
#include <algorithm>
#include <map>
#include <limits>

vtkStandardNewMacro(vtkVofSeedPoints);

namespace {
// from vtkVofAdvect
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
// from vtkVofAdvect
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
}

//-----------------------------------------------------------------------------
vtkVofSeedPoints::vtkVofSeedPoints()
{
  this->OutputSeeds = NULL;
  this->Connectivity = NULL;
  this->Coords = NULL;
}

//-----------------------------------------------------------------------------
vtkVofSeedPoints::~vtkVofSeedPoints()
{
  if (this->OutputSeeds != NULL)
    this->OutputSeeds->Delete();
  if (this->Connectivity != NULL)
    this->Connectivity->Delete();
  if (this->Coords != NULL)
    this->Coords->Delete();
}

//----------------------------------------------------------------------------
void vtkVofSeedPoints::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}

//----------------------------------------------------------------------------
void vtkVofSeedPoints::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
int vtkVofSeedPoints::FillInputPortInformation( int port, vtkInformation* info )
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

namespace {

  typedef struct {
    int x;
    int y;
    int z;
  } int3_t;

  bool int3_t_compare(const int3_t &a, const int3_t &b)
  {
    return (a.x < b.x ||
	    (a.x == b.x &&
	     (a.y < b.y ||
	      (a.y == b.y &&
	       (a.z < b.z)))));
  }

  // TODO: nodePos computation is unsafe now!
  void computeCellNodesFromCellCenters(vtkDataArray *coords,
				       vtkDataArray *coordsDer, 
				       int numCells)
  {
    float centerDist = coords->GetComponent(1,0) - coords->GetComponent(0,0);
    float nodePos = coords->GetComponent(0,0) - centerDist/2.0f;//0.0f;
    for (int i = 0; i < numCells; ++i) {
      coordsDer->SetComponent(i, 0, nodePos);
      float halfCell = coords->GetComponent(0,i) - nodePos;
      nodePos += halfCell*2.0f;
    }
    coordsDer->SetComponent(numCells, 0, nodePos); // last node
  }

  bool pointWithinBounds(const float point[3], const double bounds[6])
  {
    if (point[0] >= bounds[0] && point[0] < bounds[1] &&
	point[1] >= bounds[2] && point[1] < bounds[3] &&
	point[2] >= bounds[4] && point[2] < bounds[5]) {
      return true;
    }
    return false;
  }

  int g_seedIdx = 0;

  void placeSeeds(vtkPoints *seeds, const float cellCenter[3], 
		  const float cellSize[3], const int refinement,
		  const float f, const double bounds[6],
		  const int cell_x, const int cell_y, const int cell_z, 
		  std::map<int3_t, int, bool(*)(const int3_t &a, const int3_t &b)> &seedPos)
  {
    float originOffset[3] = {0.0f,0.0f,0.0f};
    float cellSizeTmp[3] = {cellSize[0], cellSize[1], cellSize[2]};
    int subdiv = 1;
    for (int i = 0; i < refinement; ++i) {
      cellSizeTmp[0] /= 2.0f;
      cellSizeTmp[1] /= 2.0f;
      cellSizeTmp[2] /= 2.0f;
      originOffset[0] -= cellSizeTmp[0]/2.0f;
      originOffset[1] -= cellSizeTmp[1]/2.0f;
      originOffset[2] -= cellSizeTmp[2]/2.0f;
      subdiv *= 2;
    }

    for (int zr = 0; zr < subdiv; ++zr) {
      for (int yr = 0; yr < subdiv; ++yr) {
	for (int xr = 0; xr < subdiv; ++xr) {

	  float dx[3] = {originOffset[0] + xr*cellSizeTmp[0],
			 originOffset[1] + yr*cellSizeTmp[1],
			 originOffset[2] + zr*cellSizeTmp[2]};
	  float seed[3] = {cellCenter[0]+dx[0],
			   cellCenter[1]+dx[1],
			   cellCenter[2]+dx[2]};
	  // float df = grad[0]*dx[0] + grad[1]*dx[1] + grad[2]*dx[2];

	  if (pointWithinBounds(seed, bounds)) {

	    seeds->InsertNextPoint(seed);
	    int3_t pos = {cell_x*subdiv + xr, 
			  cell_y*subdiv + yr, 
			  cell_z*subdiv + zr};
	    seedPos[pos] = g_seedIdx;
	    ++g_seedIdx;

	  }
	}
      }
    }
  }

  void computeGradient(vtkDataArray *data, const int res[3], 
		       int i, int j, int k, 
		       vtkDataArray *coordCenters[3], float grad[3])
  {
    int im = std::max(i-1,0);
    int ip = std::min(i+1,res[0]-1);
    float di = coordCenters[0]->GetComponent(ip,0) - coordCenters[0]->GetComponent(im,0);
    int jm = std::max(j-1,0);	  
    int jp = std::min(j+1,res[1]-1);
    float dj = coordCenters[1]->GetComponent(jp,0) - coordCenters[1]->GetComponent(jm,0);
    int km = std::max(k-1,0);	  
    int kp = std::min(k+1,res[2]-1);
    float dk = coordCenters[2]->GetComponent(kp,0) - coordCenters[2]->GetComponent(km,0);

    int id_left = im + j*res[0] + k*res[0]*res[1];
    int id_right = ip + j*res[0] + k*res[0]*res[1];
    int id_bottom = i + jm*res[0] + k*res[0]*res[1];
    int id_top = i + jp*res[0] + k*res[0]*res[1];
    int id_back = i + j*res[0] + km*res[0]*res[1];
    int id_front = i + j*res[0] + kp*res[0]*res[1];

    grad[0] = (data->GetComponent(id_right,0) - 
	       data->GetComponent(id_left,0))/di;
    grad[1] = (data->GetComponent(id_top,0) - 
	       data->GetComponent(id_bottom,0))/dj;
    grad[2] = (data->GetComponent(id_front,0) - 
	       data->GetComponent(id_back,0))/dk;
  }
} // namespace
int vtkVofSeedPoints::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
					  vtkInformationVector **inputVector,
					  vtkInformationVector *outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);

  if (Reseed) {
    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {

      unsigned int numberOfInputTimeSteps =
	inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

      std::vector<double> inputTimeValues(numberOfInputTimeSteps);
      inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
		  &inputTimeValues[0]);

      inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), 
		  inputTimeValues[SeedTimeStep]);
    }
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkVofSeedPoints::RequestData(vtkInformation *request,
				  vtkInformationVector **inputVector,
				  vtkInformationVector *outputVector)
{
  if (!(this->OutputSeeds == NULL || Reseed)) {

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    output->SetPoints(OutputSeeds);
    output->GetPointData()->AddArray(Connectivity);
    output->GetPointData()->AddArray(Coords);

    return 1;
  }

  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  // get the input
  vtkRectilinearGrid *input = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // determine if data is cell-based or point-based
  vtkDataArray *pointData = input->GetPointData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  vtkDataArray *cellData = input->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  vtkDataArray *data;

  int inputRes[3];
  input->GetDimensions(inputRes);
  int cellRes[3] = {inputRes[0], inputRes[1], inputRes[2]};

  vtkDataArray *coordCenters[3];
  vtkDataArray *coordNodes[3];

  if (pointData == 0 && cellData != 0) {  // cell data
    data = cellData;
    cellRes[0] -= 1;
    cellRes[1] -= 1;
    cellRes[2] -= 1;
    DataOnCells = true;

    coordNodes[0] = input->GetXCoordinates();
    coordNodes[1] = input->GetYCoordinates();
    coordNodes[2] = input->GetZCoordinates();

    for (int c = 0; c < 3; ++c) {
      coordCenters[c] = vtkFloatArray::New();
      coordCenters[c]->SetNumberOfComponents(1);
      coordCenters[c]->SetNumberOfTuples(coordNodes[c]->GetNumberOfTuples()-1);
      for (int i = 0; i < coordCenters[c]->GetNumberOfTuples(); ++i) {
	coordCenters[c]->SetComponent(i,0,(coordNodes[c]->GetComponent(0,i) + 
					   coordNodes[c]->GetComponent(0,i+1))/2.0f);
      }
    }
  }
  else if (pointData != 0 && cellData == 0) {  // point data
    data = pointData;
    DataOnCells = false;

    coordCenters[0] = input->GetXCoordinates();
    coordCenters[1] = input->GetYCoordinates();
    coordCenters[2] = input->GetZCoordinates();

    for (int c = 0; c < 3; ++c) {
      coordNodes[c] = vtkFloatArray::New();
      coordNodes[c]->SetNumberOfComponents(1);
      coordNodes[c]->SetNumberOfTuples(cellRes[c]+1);
      computeCellNodesFromCellCenters(coordCenters[c], coordNodes[c], cellRes[c]);
    }
  }

  if (OutputSeeds != NULL) {
    OutputSeeds->Delete();
  }
  OutputSeeds = vtkPoints::New();
  OutputSeeds->SetDataTypeToFloat();

  double bounds[6];
  input->GetBounds(bounds);
  int extent[6];
  input->GetExtent(extent);

  std::map<int3_t, int, bool(*)(const int3_t &a, const int3_t &b)> seedPos(int3_t_compare);
  seedPos.clear();
  g_seedIdx = 0;

  //---------------------------------------------------------------------------
  // populate the grid with seed points
  int idx = 0;
  for (int k = 0; k < cellRes[2]; ++k) {
    for (int j = 0; j < cellRes[1]; ++j) {
      for (int i = 0; i < cellRes[0]; ++i) {

	float f = data->GetComponent(0,idx);
	if (f > 0.0f) {
	  float cellCenter[3] = {coordCenters[0]->GetComponent(i,0),
				 coordCenters[1]->GetComponent(j,0),
				 coordCenters[2]->GetComponent(k,0)};
	  float cellSize[3] = {coordNodes[0]->GetComponent(i+1,0) - 
			       coordNodes[0]->GetComponent(i,0),
			       coordNodes[1]->GetComponent(j+1,0) - 
			       coordNodes[1]->GetComponent(j,0),
			       coordNodes[2]->GetComponent(k+1,0) - 
			       coordNodes[2]->GetComponent(k,0)};

	  placeSeeds(OutputSeeds, cellCenter, cellSize, Refinement, f, 
		     bounds, i+extent[0], j+extent[2], k+extent[4], seedPos);
	}
	++idx;
      }
    }
  }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  output->SetPoints(OutputSeeds);

  if (Connectivity != NULL) {
    Connectivity->Delete();
  }
  Connectivity = vtkIntArray::New();
  Connectivity->SetName("Connectivity");
  Connectivity->SetNumberOfComponents(3);
  Connectivity->SetNumberOfTuples(seedPos.size());

  if (Coords != NULL) {
    Coords->Delete();
  }
  Coords = vtkShortArray::New();
  Coords->SetName("Coords");
  Coords->SetNumberOfComponents(3);
  Coords->SetNumberOfTuples(seedPos.size());

  std::map<int3_t, int>::iterator it;
  for (it = seedPos.begin(); it != seedPos.end(); ++it) {

    int conn[3] = {-1,-1,-1};
    int x = it->first.x;
    int y = it->first.y;
    int z = it->first.z;

    int3_t pos_xm = {x-1,y,z};
    if (seedPos.find(pos_xm) != seedPos.end()) {
      conn[0] = seedPos[pos_xm];
    }
    int3_t pos_ym = {x,y-1,z};
    if (seedPos.find(pos_ym) != seedPos.end()) {
      conn[1] = seedPos[pos_ym];
    }
    int3_t pos_zm = {x,y,z-1};
    if (seedPos.find(pos_zm) != seedPos.end()) {
      conn[2] = seedPos[pos_zm];
    }
    
    int seedIdx = it->second;
    Connectivity->SetTuple3(seedIdx, conn[0], conn[1], conn[2]);
    Coords->SetTuple3(seedIdx, x, y, z);
  }

  if (Connectivity != NULL) {
    output->GetPointData()->AddArray(Connectivity);
  }
  if (Coords != NULL) {
    output->GetPointData()->AddArray(Coords);
  }

  // // used only for testing ---------------------------------------------------
  // vtkCellArray *lines = vtkCellArray::New();
  // for (int i = 0; i < Connectivity->GetNumberOfTuples(); ++i) {
  //   int conn[3] = {Connectivity->GetComponent(i,0),
  // 		   Connectivity->GetComponent(i,1),
  // 		   Connectivity->GetComponent(i,2)};
  //   vtkIdType pts[2];
  //   pts[0] = i;

  //   pts[1] = conn[0];
  //   if (pts[1] != -1) {
  //     lines->InsertNextCell(2, pts);
  //   }
  //   pts[1] = conn[1];
  //   if (pts[1] != -1) {
  //     lines->InsertNextCell(2, pts);
  //   }
  //   pts[1] = conn[2];
  //   if (pts[1] != -1) {
  //     lines->InsertNextCell(2, pts);
  //   }
  // }
  // output->SetLines(lines);

  return 1;
}

////////// External Operators /////////////
void vtkVofSeedPoints::PrintSelf(ostream &os, vtkIndent indent)
{
}
