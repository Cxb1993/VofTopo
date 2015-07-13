#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkShortArray.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkVofGenBounds.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <limits>

vtkStandardNewMacro(vtkVofGenBounds);

namespace {
  typedef struct {
    int x;
    int y;
    int z;
  } int3_t;

  typedef struct {
    float x;
    float y;
    float z;
  } float3_t;

  void operator+=(float3_t &a, float3_t b)
  {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
  }

  float3_t operator/(float3_t a, float b)
  {
    float3_t c = {a.x/b, a.y/b, a.z/b};
    return c;
  }

  float3_t operator-(float3_t a, float3_t b)
  {
    float3_t c = {a.x-b.x, a.y-b.y, a.z-b.z};
    return c;
  }

  float3_t cross(float3_t a, float3_t b)
  {
    float3_t c = {a.y*b.z - a.z*b.y, 
		  a.z*b.x - a.x*b.z, 
		  a.x*b.y - a.y*b.x};
    return c;
  }

  float length(float3_t a)
  {
    return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  }

  float3_t normalize(float3_t a)
  {
    float len = length(a);
    if (len > 0.0f) 
      return a/len;
    float3_t z = {0.0f,0.0f,0.0f};
    return z;
  }

  class compare_int3_t {
  public:
    bool operator()(const int3_t a, const int3_t b) const {
      return (a.x < b.x || (a.x == b.x && (a.y < b.y || (a.y == b.y && (a.z < b.z)))));
    }
  };

  void mergeTriangles(std::vector<float3_t>& vertices,
		      std::vector<int3_t>& ivertices,
		      std::vector<int>& indices,
		      std::vector<float3_t>& mergedVertices)
  {
    int vertexID = 0;
    std::map<int3_t, int, compare_int3_t> vertexMap;
    int totalVerts = vertices.size();
    
    for (int t = 0; t < totalVerts; t++) {

      int3_t &key = ivertices[t];
      
      if (vertexMap.find(key) == vertexMap.end()) {
	
	vertexMap[key] = vertexID;

	mergedVertices.push_back(vertices[t]);

	indices.push_back(vertexID);
	
	vertexID++;
      }
      else {	
	indices.push_back(vertexMap[key]);
      }
    }
  }

  void smoothSurface(std::vector<float3_t>& vertices,
		     std::vector<int>& indices)
  {
    std::vector<std::set<int> > neighbors(vertices.size());
    for (int i = 0; i < indices.size()/3; ++i) {
      
      int *tri = &indices[i*3];
      for (int j = 0; j < 3; j++) {
	neighbors[tri[j]].insert(tri[(j+1)%3]);
	neighbors[tri[j]].insert(tri[(j+2)%3]);
      }
    }

    std::vector<float3_t> verticesTmp(vertices.size());
    float3_t vert = {0.0f,0.0f,0.0f};
    for (int i = 0; i < verticesTmp.size(); ++i) {
      verticesTmp[i] = vert;
    }

    for (int i = 0; i < neighbors.size(); ++i) {
      std::set<int>::iterator it;
      for (it = neighbors[i].begin(); it != neighbors[i].end(); ++it) {
	verticesTmp[i] += vertices[*it];
      }
      verticesTmp[i] += vertices[i];
      vertices[i] = verticesTmp[i]/(neighbors[i].size()+1);
    }
  }

  void generateNormals(std::vector<float3_t>& vertices,
		       std::vector<int>& indices,
		       std::vector<float3_t>& normals)
  {
    normals.resize(vertices.size());
    float3_t norm = {0.0f,0.0f,0.0f};
    for (int i = 0; i < normals.size(); ++i) {
      normals[i] = norm;
    }

    for (int i = 0; i < indices.size()/3; ++i) {
      
      int *tri = &indices[i*3];
      float3_t v0 = vertices[tri[0]];
      float3_t v1 = vertices[tri[1]];
      float3_t v2 = vertices[tri[2]];

      float3_t e0 = v1 - v0;
      float3_t e1 = v2 - v0;
      float3_t cr = cross(e0,e1);
      
      normals[tri[0]] += cr;
      normals[tri[1]] += cr;
      normals[tri[2]] += cr;
    }
    for (int i = 0; i < normals.size(); ++i) {
      normals[i] = normalize(normals[i]);
    }
  }
}

//----------------------------------------------------------------------------
vtkVofGenBounds::vtkVofGenBounds()
{
}

//----------------------------------------------------------------------------
vtkVofGenBounds::~vtkVofGenBounds()
{
}

//----------------------------------------------------------------------------
void vtkVofGenBounds::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}

//----------------------------------------------------------------------------
void vtkVofGenBounds::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
int vtkVofGenBounds::FillInputPortInformation( int port, vtkInformation* info )
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
int vtkVofGenBounds::RequestData(vtkInformation *request,
				 vtkInformationVector **inputVector,
				 vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPoints *points = input->GetPoints();
  vtkFloatArray *labels = vtkFloatArray::
    SafeDownCast(input->GetPointData()->GetArray("Labels"));
  vtkIntArray *connectivity = vtkIntArray::
    SafeDownCast(input->GetPointData()->GetArray("Connectivity"));
  vtkShortArray *coords = vtkShortArray::
    SafeDownCast(input->GetPointData()->GetArray("Coords"));

  const int numPoints = points->GetNumberOfPoints();

  // hack: compute cell size - should actually be taken from the grid
  // find a seed point that has neighbors in x, y, and z directions - from 
  // distance to these points we can compute the cell size; without the hack
  // it must be taken from grid explicitely
  float cellSize[3];

  for (int i = 0; i < numPoints; ++i) {

    int conn[3] = {connectivity->GetComponent(i, 0),
		   connectivity->GetComponent(i, 1),
		   connectivity->GetComponent(i, 2)};
    if (conn[0] > -1 && conn[1] > -1 && conn[2] > -1) {
      
      double p0[3];
      points->GetPoint(i, p0);
      double p1[3];
      points->GetPoint(conn[0], p1);
      double p2[3];
      points->GetPoint(conn[1], p2);
      double p3[3];
      points->GetPoint(conn[2], p3);

      cellSize[0] = std::abs(p1[0]-p0[0]);
      cellSize[1] = std::abs(p2[1]-p0[1]);
      cellSize[2] = std::abs(p3[2]-p0[2]);
      break;
    }
  }
  // hack end
  const float cs2[3] = {cellSize[0]/2.0f, cellSize[1]/2.0f, cellSize[2]/2.0f};

  std::vector<int3_t> ivertices;
  ivertices.clear();
  std::vector<float3_t> vertices;
  vertices.clear();

  const int co[6][12] = {{0,0,0, 0,0,1, 0,1,1, 0,1,0},
  			 {0,0,0, 1,0,0, 1,0,1, 0,0,1},
  			 {0,0,0, 0,1,0, 1,1,0, 1,0,0},
			 {1,0,0, 1,0,1, 1,1,1, 1,1,0},
  			 {0,1,0, 1,1,0, 1,1,1, 0,1,1},
  			 {0,0,1, 0,1,1, 1,1,1, 1,0,1}};
  const float po[6][12] = {{-cs2[0],-cs2[1],-cs2[2], 
  			    -cs2[0],-cs2[1], cs2[2], 
  			    -cs2[0], cs2[1], cs2[2], 
  			    -cs2[0], cs2[1],-cs2[2]},
  			   {-cs2[0],-cs2[1],-cs2[2], 
  			     cs2[0],-cs2[1],-cs2[2], 
  			     cs2[0],-cs2[1], cs2[2], 
  			    -cs2[0],-cs2[1], cs2[2]},
  			   {-cs2[0],-cs2[1],-cs2[2], 
  			    -cs2[0], cs2[1],-cs2[2], 
  			     cs2[0], cs2[1],-cs2[2],
			     cs2[0],-cs2[1],-cs2[2]},
			   { cs2[0],-cs2[1],-cs2[2], 
  			     cs2[0],-cs2[1], cs2[2], 
  			     cs2[0], cs2[1], cs2[2], 
  			     cs2[0], cs2[1],-cs2[2]},
  			   {-cs2[0], cs2[1],-cs2[2], 
  			     cs2[0], cs2[1],-cs2[2], 
  			     cs2[0], cs2[1], cs2[2], 
  			    -cs2[0], cs2[1], cs2[2]},
  			   {-cs2[0],-cs2[1], cs2[2], 
  			    -cs2[0], cs2[1], cs2[2], 
  			     cs2[0], cs2[1], cs2[2],
  			     cs2[0],-cs2[1], cs2[2]}};

  vtkFloatArray *labelsBack = vtkFloatArray::New();
  labelsBack->SetNumberOfComponents(1);
  labelsBack->SetName("BackLabels");
  vtkFloatArray *labelsFront = vtkFloatArray::New();
  labelsFront->SetNumberOfComponents(1);
  labelsFront->SetName("FrontLabels");

  // std::vector<std::vector<int> > neighbors(numPoints);
  // std::vector<int> initNeighbors(6,-1);
  // for (int i = 0; i < numPoints; ++i) {
  //   neighbors[i] = initNeighbors;
  // }

  for (int i = 0; i < numPoints; ++i) {

    int conn[3] = {connectivity->GetComponent(i, 0),
		   connectivity->GetComponent(i, 1),
		   connectivity->GetComponent(i, 2)};

    float l0 = labels->GetValue(i);
    int c0[3] = {coords->GetComponent(i, 0),
		 coords->GetComponent(i, 1),
		 coords->GetComponent(i, 2)};

    double p0[3];
    points->GetPoint(i, p0);

    for (int j = 0; j < 3; ++j) {
      if (conn[j] > -1) {

	// neighbors[i][j] = 1;
	// neighbors[conn[j]][j+3] = 1;

	float l1 = labels->GetValue(conn[j]);
	double p1[3];
	points->GetPoint(conn[j], p1);
      
	if (l0 != l1) {

	  labelsBack->InsertNextValue(l0);
	  labelsBack->InsertNextValue(l0);
	  labelsFront->InsertNextValue(l1);
	  labelsFront->InsertNextValue(l1);

	  float3_t verts[4];
	  int3_t iverts[4];

	  for (int k = 0; k < 4; ++k) {

	    float3_t vertex = {p0[0]+po[j][k*3+0], 
			       p0[1]+po[j][k*3+1], 
			       p0[2]+po[j][k*3+2]};
	    verts[k] = vertex;

	    int3_t ivertex = {c0[0]+co[j][k*3+0], 
	    		      c0[1]+co[j][k*3+1], 
	    		      c0[2]+co[j][k*3+2]};
	    iverts[k] = ivertex;
	  }
	  vertices.push_back(verts[0]);
	  vertices.push_back(verts[1]);
	  vertices.push_back(verts[2]);
	  vertices.push_back(verts[2]);
	  vertices.push_back(verts[3]);
	  vertices.push_back(verts[0]);

	  ivertices.push_back(iverts[0]);
	  ivertices.push_back(iverts[1]);
	  ivertices.push_back(iverts[2]);
	  ivertices.push_back(iverts[2]);
	  ivertices.push_back(iverts[3]);
	  ivertices.push_back(iverts[0]);
	}
      }
    }    
  }

  // float minval = -1.0f;//-1.0f*std::numeric_limits<float>::max();

  // for (int i = 0; i < numPoints; ++i) {

  //   double p0[3];
  //   points->GetPoint(i, p0);
  //   int c0[3] = {coords->GetComponent(i, 0),
  // 		 coords->GetComponent(i, 1),
  // 		 coords->GetComponent(i, 2)};
  //   float l0 = labels->GetValue(i);

  //   for (int j = 0; j < 6; ++j) {
  //     if (neighbors[i][j] == -1) {


  // 	if (j < 3) {
  // 	  labelsBack->InsertNextValue(l0);
  // 	  labelsBack->InsertNextValue(l0);
  // 	  labelsFront->InsertNextValue(minval);
  // 	  labelsFront->InsertNextValue(minval);
  // 	}
  // 	else {
  // 	  labelsBack->InsertNextValue(minval);
  // 	  labelsBack->InsertNextValue(minval);
  // 	  labelsFront->InsertNextValue(l0);
  // 	  labelsFront->InsertNextValue(l0);
  // 	}
	
  // 	float3_t verts[4];
  // 	int3_t iverts[4];

  // 	for (int k = 0; k < 4; ++k) {

  // 	  float3_t vertex = {p0[0]+po[j][k*3+0], 
  // 			     p0[1]+po[j][k*3+1], 
  // 			     p0[2]+po[j][k*3+2]};
  // 	  verts[k] = vertex;

  // 	  int3_t ivertex = {c0[0]+co[j][k*3+0], 
  // 			    c0[1]+co[j][k*3+1], 
  // 			    c0[2]+co[j][k*3+2]};
  // 	  iverts[k] = ivertex;
  // 	}
  // 	vertices.push_back(verts[0]);
  // 	vertices.push_back(verts[1]);
  // 	vertices.push_back(verts[2]);
  // 	vertices.push_back(verts[2]);
  // 	vertices.push_back(verts[3]);
  // 	vertices.push_back(verts[0]);

  // 	ivertices.push_back(iverts[0]);
  // 	ivertices.push_back(iverts[1]);
  // 	ivertices.push_back(iverts[2]);
  // 	ivertices.push_back(iverts[2]);
  // 	ivertices.push_back(iverts[3]);
  // 	ivertices.push_back(iverts[0]);
  //     }
  //   }
  // }

  std::vector<int> indices;
  std::vector<float3_t> mergedVertices;
  mergeTriangles(vertices, ivertices, indices, mergedVertices);
  smoothSurface(mergedVertices, indices);
  
  vtkPoints *outputPoints = vtkPoints::New();
  outputPoints->SetNumberOfPoints(mergedVertices.size());
  for (int i = 0; i < mergedVertices.size(); ++i) {
    
    double p[3] = {mergedVertices[i].x,
		   mergedVertices[i].y,
		   mergedVertices[i].z};
    outputPoints->SetPoint(i, p);        
  }

  vtkIdTypeArray *cells = vtkIdTypeArray::New();
  cells->SetNumberOfComponents(1);
  cells->SetNumberOfTuples(indices.size()/3*4);
  for (int i = 0; i < indices.size()/3; ++i) {
    cells->SetValue(i*4+0,3);
    cells->SetValue(i*4+1,indices[i*3+0]);
    cells->SetValue(i*4+2,indices[i*3+1]);
    cells->SetValue(i*4+3,indices[i*3+2]);
  }

  vtkCellArray *outputTriangles = vtkCellArray::New();
  outputTriangles->SetNumberOfCells(indices.size()/3);
  outputTriangles->SetCells(indices.size()/3, cells);

  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  output->SetPoints(outputPoints);
  output->SetPolys(outputTriangles);

  output->GetCellData()->AddArray(labelsBack);
  output->GetCellData()->AddArray(labelsFront);

  return 1;
}

////////// External Operators /////////////
void vtkVofGenBounds::PrintSelf(ostream &os, vtkIndent indent)
{
}
