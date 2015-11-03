#include "vofTopology.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <set>
#include <cmath>

namespace
{

  bool compare_int3(const int3 &a, const int3 &b)
  {
    return (a.x < b.x ||
  	    (a.x == b.x &&
  	     (a.y < b.y ||
  	      (a.y == b.y &&
  	       (a.z < b.z)))));
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

  void placeSeeds(vtkPoints *seeds, const float cellCenter[3], 
		  const float cellSize[3], const int refinement,
		  const float f, const float gradf[3], const double bounds[6],
		  const int cell_x, const int cell_y, const int cell_z, 
		  std::map<int3, int, bool(*)(const int3 &a, const int3 &b)> &seedPos,
		  int &seedIdx)
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
	  float df = gradf[0]*dx[0] + gradf[1]*dx[1] + gradf[2]*dx[2];

	  if (pointWithinBounds(seed, bounds) && f + df >= 0.06125f) {
	    // if (pointWithinBounds(seed, bounds) && f >= 0.99f) {

	    seeds->InsertNextPoint(seed);
	    int3 pos = {cell_x*subdiv + xr, 
			  cell_y*subdiv + yr, 
			  cell_z*subdiv + zr};
	    seedPos[pos] = seedIdx;
	    ++seedIdx;

	  }
	}
      }
    }    
  }

  void placeSeedsPLIC(vtkPoints *seeds,
		      const float cellCenter[3], 
		      const float cellSize[3],
		      const int refinement,
		      const int cellRes[3],
		      const float f,
		      const std::vector<float> &lstar,
		      const std::vector<float> &normalsInt,
		      const double bounds[6],
		      const int cell_x, const int cell_y, const int cell_z,
		      const int idx,
		      std::map<int3, int, bool(*)(const int3 &a, const int3 &b)> &seedPos,
		      int &seedIdx)
  {
    float originOffset[3] = {0.0f,0.0f,0.0f};
    float cellSizeTmp[3] = {cellSize[0], cellSize[1], cellSize[2]};
    int subdiv = 1;
    float offset[3] = {0.5f,0.5f,0.5f};
    float scale[3] = {1.0f,1.0f,1.0f}; 
    for (int i = 0; i < refinement; ++i) {
      cellSizeTmp[0] /= 2.0f;
      cellSizeTmp[1] /= 2.0f;
      cellSizeTmp[2] /= 2.0f;
      originOffset[0] -= cellSizeTmp[0]/2.0f;
      originOffset[1] -= cellSizeTmp[1]/2.0f;
      originOffset[2] -= cellSizeTmp[2]/2.0f;
      subdiv *= 2;

      offset[0] /= 2.0f;
      offset[1] /= 2.0f;
      offset[2] /= 2.0f;
      scale[0] /= 2.0f;
      scale[1] /= 2.0f;
      scale[2] /= 2.0f;
    }

    float attachPoint[3] =
      {normalsInt[idx*3+0]>0 ? cellCenter[0]-cellSize[0]/2.0f : cellCenter[0]+cellSize[0]/2.0f,
       normalsInt[idx*3+1]>0 ? cellCenter[1]-cellSize[1]/2.0f : cellCenter[1]+cellSize[1]/2.0f,
       normalsInt[idx*3+2]>0 ? cellCenter[2]-cellSize[2]/2.0f : cellCenter[2]+cellSize[2]/2.0f};
    float n[3] = {normalsInt[idx*3+0],
		  normalsInt[idx*3+1],
		  normalsInt[idx*3+2]};

    for (int zr = 0; zr < subdiv; ++zr) {
      for (int yr = 0; yr < subdiv; ++yr) {
	for (int xr = 0; xr < subdiv; ++xr) {

	  float dx[3] = {originOffset[0] + xr*cellSizeTmp[0],
			 originOffset[1] + yr*cellSizeTmp[1],
			 originOffset[2] + zr*cellSizeTmp[2]};
	  float seed[3] = {cellCenter[0]+dx[0], // pos within cell
			   cellCenter[1]+dx[1],
			   cellCenter[2]+dx[2]};

	  float posVec[3] = {seed[0] - attachPoint[0],
			     seed[1] - attachPoint[1],
			     seed[2] - attachPoint[2]};
	  float d = posVec[0]*n[0] + posVec[1]*n[1] + posVec[2]*n[2];
	  d = std::abs(d);

	  if (pointWithinBounds(seed, bounds) &&
	      (f < g_emf1 && d < lstar[idx] || f >= g_emf1)) {	    

	    seeds->InsertNextPoint(seed);
	    int3 pos = {cell_x*subdiv + xr, 
			cell_y*subdiv + yr, 
			cell_z*subdiv + zr};
	    seedPos[pos] = seedIdx;
	    ++seedIdx;
	  }
	}
      }
    }    
  }

  static float interpolateSca(vtkDataArray *vofField,
			      const int* res, const int idxCell[3],
			      const double bcoords[3])
  {
    int lx = idxCell[0];
    int ly = idxCell[1];
    int lz = idxCell[2];
    float x = bcoords[0] - 0.5;
    float y = bcoords[1] - 0.5;
    float z = bcoords[2] - 0.5;

    if (bcoords[0] < 0.5) {
      lx -= 1;
      x = bcoords[0] + 0.5;
    }
    if (bcoords[1] < 0.5) {
      ly -= 1;
      y = bcoords[1] + 0.5;
    }
    if (bcoords[2] < 0.5) {
      lz -= 1;
      z = bcoords[2] + 0.5;
    }
    
    int ux = lx+1;
    int uy = ly+1;
    int uz = lz+1;

    if (lx < 0) lx = 0;
    if (ly < 0) ly = 0;
    if (lz < 0) lz = 0;
    if (ux > res[0]-1) ux = res[0]-1;
    if (uy > res[1]-1) uy = res[1]-1;
    if (uz > res[2]-1) uz = res[2]-1;

    unsigned lzslab = lz*res[0]*res[1];
    unsigned uzslab = uz*res[0]*res[1];
    int lyr = ly*res[0];
    int uyr = uy*res[0];

    unsigned id[8] = {lx + lyr + lzslab,
		      ux + lyr + lzslab,
		      lx + uyr + lzslab,
		      ux + uyr + lzslab,
		      lx + lyr + uzslab,
		      ux + lyr + uzslab,
		      lx + uyr + uzslab,
		      ux + uyr + uzslab};
    float vv[8];
    for (int i = 0; i < 8; i++) {
      vv[i] = vofField->GetComponent(id[i], 0);
    }

    float a = (1.0f-x)*vv[0] + x*vv[1];
    float b = (1.0f-x)*vv[2] + x*vv[3];
    float c = (1.0f-y)*a + y*b;
    a = (1.0f-x)*vv[4] + x*vv[5];
    b = (1.0f-x)*vv[6] + x*vv[7];
    float d = (1.0f-y)*a + y*b;

    return (1.0f-z)*c + z*d;
  }

  void computeGradient(vtkRectilinearGrid *grid, vtkDataArray *data,
		       const int res[3], int ijk[3],
		       vtkDataArray *coordCenters[3],
		       double pcoords[3], float grad[3])
  {
    int i = ijk[0];
    int j = ijk[1];
    int k = ijk[2];
    int im = std::max(i-1,0);
    int ip = std::min(i+1,res[0]-1);
    float dx = coordCenters[0]->GetComponent(ip,0) - coordCenters[0]->GetComponent(im,0);
    int jm = std::max(j-1,0);	  
    int jp = std::min(j+1,res[1]-1);
    float dy = coordCenters[1]->GetComponent(jp,0) - coordCenters[1]->GetComponent(jm,0);
    int km = std::max(k-1,0);	  
    int kp = std::min(k+1,res[2]-1);
    float dz = coordCenters[2]->GetComponent(kp,0) - coordCenters[2]->GetComponent(km,0);

    int id_left   = im + j*res[0] + k*res[0]*res[1];
    int id_right  = ip + j*res[0] + k*res[0]*res[1];
    int id_bottom = i + jm*res[0] + k*res[0]*res[1];
    int id_top    = i + jp*res[0] + k*res[0]*res[1];
    int id_back   = i + j*res[0] + km*res[0]*res[1];
    int id_front  = i + j*res[0] + kp*res[0]*res[1];

    // int f_left = ;
    // int f_right;
    // int f_bottom;
    // int f_top;   
    // int f_back;
    // int f_front;
    
    grad[0] = (data->GetComponent(id_right,0) - 
	       data->GetComponent(id_left,0))/dx;
    grad[1] = (data->GetComponent(id_top,0) - 
	       data->GetComponent(id_bottom,0))/dy;
    grad[2] = (data->GetComponent(id_front,0) - 
	       data->GetComponent(id_back,0))/dz;
  }

  void computeGradient(vtkDataArray *data, const int res[3], 
		       int i, int j, int k, 
		       vtkDataArray *coordCenters[3], float grad[3])
  {
    int im = std::max(i-1,0);
    int ip = std::min(i+1,res[0]-1);
    float dx = coordCenters[0]->GetComponent(ip,0) - coordCenters[0]->GetComponent(im,0);
    int jm = std::max(j-1,0);	  
    int jp = std::min(j+1,res[1]-1);
    float dy = coordCenters[1]->GetComponent(jp,0) - coordCenters[1]->GetComponent(jm,0);
    int km = std::max(k-1,0);	  
    int kp = std::min(k+1,res[2]-1);
    float dz = coordCenters[2]->GetComponent(kp,0) - coordCenters[2]->GetComponent(km,0);

    int id_left = im + j*res[0] + k*res[0]*res[1];
    int id_right = ip + j*res[0] + k*res[0]*res[1];
    int id_bottom = i + jm*res[0] + k*res[0]*res[1];
    int id_top = i + jp*res[0] + k*res[0]*res[1];
    int id_back = i + j*res[0] + km*res[0]*res[1];
    int id_front = i + j*res[0] + kp*res[0]*res[1];

    grad[0] = (data->GetComponent(id_right,0) - 
	       data->GetComponent(id_left,0))/dx;
    grad[1] = (data->GetComponent(id_top,0) - 
	       data->GetComponent(id_bottom,0))/dy;
    grad[2] = (data->GetComponent(id_front,0) - 
	       data->GetComponent(id_back,0))/dz;
  }

  static float3 interpolateVec(vtkDataArray *velocityField,
			       const int* res, const int idxCell[3],
			       const double bcoords[3])
  {
    int lx = idxCell[0];
    int ly = idxCell[1];
    int lz = idxCell[2];
    float x = bcoords[0] - 0.5;
    float y = bcoords[1] - 0.5;
    float z = bcoords[2] - 0.5;

    if (bcoords[0] < 0.5) {
      lx -= 1;
      x = bcoords[0] + 0.5;
    }
    if (bcoords[1] < 0.5) {
      ly -= 1;
      y = bcoords[1] + 0.5;
    }
    if (bcoords[2] < 0.5) {
      lz -= 1;
      z = bcoords[2] + 0.5;
    }
    
    int ux = lx+1;
    int uy = ly+1;
    int uz = lz+1;

    if (lx < 0) lx = 0;
    if (ly < 0) ly = 0;
    if (lz < 0) lz = 0;
    if (ux > res[0]-1) ux = res[0]-1;
    if (uy > res[1]-1) uy = res[1]-1;
    if (uz > res[2]-1) uz = res[2]-1;

    unsigned lzslab = lz*res[0]*res[1];
    unsigned uzslab = uz*res[0]*res[1];
    int lyr = ly*res[0];
    int uyr = uy*res[0];

    unsigned id[8] = {lx + lyr + lzslab,
		      ux + lyr + lzslab,
		      lx + uyr + lzslab,
		      ux + uyr + lzslab,
		      lx + lyr + uzslab,
		      ux + lyr + uzslab,
		      lx + uyr + uzslab,
		      ux + uyr + uzslab};
    float3 vv[8];
    for (int i = 0; i < 8; i++) {
      vv[i].x = velocityField->GetComponent(id[i], 0);
      vv[i].y = velocityField->GetComponent(id[i], 1);
      vv[i].z = velocityField->GetComponent(id[i], 2);
    }

    float3 a = (1.0f-x)*vv[0] + x*vv[1];
    float3 b = (1.0f-x)*vv[2] + x*vv[3];
    float3 c = (1.0f-y)*a + y*b;
    a = (1.0f-x)*vv[4] + x*vv[5];
    b = (1.0f-x)*vv[6] + x*vv[7];
    float3 d = (1.0f-y)*a + y*b;

    return (1.0f-z)*c + z*d;
  }

  // connected-components
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

  class compare_float3 {
  public:
    bool operator()(const float3 a, const float3 b) const {
      return (a.x < b.x || (a.x == b.x && (a.y < b.y || (a.y == b.y && (a.z < b.z)))));
    }
  };

  void mergeTriangles(std::vector<float3>& vertices,
		      std::vector<float3>& ivertices,
		      std::vector<int>& indices,
		      std::vector<float3>& mergedVertices)
  {
    int vertexID = 0;
    std::map<float3, int, compare_float3> vertexMap;
    int totalVerts = vertices.size();
    
    for (int t = 0; t < totalVerts; t++) {

      float3 &key = ivertices[t];
      
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

  void mergeVertices(const std::vector<float3>& ivertices,
  		     const int vertexOffset,
  		     std::vector<int>& indices,
  		     std::vector<int>& vertexList)
  {
    std::map<float3, int, compare_float3> vertexMap;
    int totalVerts = ivertices.size();
    int vertexID = vertexOffset;
    
    for (int t = 0; t < totalVerts; t++) {

      const float3 &key = ivertices[t];
      
      if (vertexMap.find(key) == vertexMap.end()) {
	
  	vertexMap[key] = vertexID;
  	vertexList.push_back(t);
  	indices.push_back(vertexID);

  	vertexID++;	
      }
      else {	
  	indices.push_back(vertexMap[key]);
      }
    }
  }
  void mergeVertices(const std::vector<float3>& ivertices,
		     const int vertexOffset,
		     vtkCellArray *cells,
		     std::vector<int>& vertexList)
  {
    std::map<float3, int, compare_float3> vertexMap;
    int totalVerts = ivertices.size();
    int numTriangles = totalVerts/3;
    int vertexID = vertexOffset;
    
    for (int t = 0; t < numTriangles; t++) {
      vtkIdType pts[3];
      for (int v = 0; v < 3; ++v) {

	const float3 &key = ivertices[t*3+v];
      
	if (vertexMap.find(key) == vertexMap.end()) {
	
	  vertexMap[key] = vertexID;
	  vertexList.push_back(t*3+v);
	  pts[v] = vertexID;
	  vertexID++;	
	}
	else {	
	  pts[v] = vertexMap[key];
	}
      }
      cells->InsertNextCell(3,pts);      
    }
  }

}

// taken from vtkParticleTracerBase.cxx
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

void smoothSurface(std::vector<float3>& vertices,
		   std::vector<int>& indices)
{
  std::vector<std::set<int> > neighbors(vertices.size());
  for (int i = 0; i < indices.size()/3; ++i) {
      
    int *tri = &indices[i*3];
    for (int j = 0; j < 3; j++) {
      neighbors[tri[j]].insert(tri[(j+0)%3]);
      neighbors[tri[j]].insert(tri[(j+1)%3]);
      neighbors[tri[j]].insert(tri[(j+2)%3]);
    }
  }

  std::vector<float3> verticesTmp(vertices.size());
  float3 vert = {0.0f,0.0f,0.0f};
  for (int i = 0; i < verticesTmp.size(); ++i) {
    verticesTmp[i] = vert;
  }

  for (int i = 0; i < neighbors.size(); ++i) {
    std::set<int>::iterator it;
    for (it = neighbors[i].begin(); it != neighbors[i].end(); ++it) {
      verticesTmp[i] += vertices[*it];
    }
  }

  for (int i = 0; i < verticesTmp.size(); ++i) {
    if (neighbors.size() > 5) {
      vertices[i] = verticesTmp[i]/(neighbors[i].size());
    }
  }
}

void smoothSurface(vtkPoints *vertices,
		   vtkCellArray *cells)
{
  int numVertices = vertices->GetNumberOfPoints();
  std::vector<std::set<int> > neighbors(numVertices);
  vtkIdTypeArray *indices = cells->GetData();
  
  for (int i = 0; i < cells->GetNumberOfCells(); ++i) {
      
    int tri[3] = {indices->GetTuple1(i*4 + 1),
		  indices->GetTuple1(i*4 + 2),
		  indices->GetTuple1(i*4 + 3)};
    for (int j = 0; j < 3; j++) {
      neighbors[tri[j]].insert(tri[(j+0)%3]);
      neighbors[tri[j]].insert(tri[(j+1)%3]);
      neighbors[tri[j]].insert(tri[(j+2)%3]);
    }
  }

  std::vector<float3> verticesTmp(numVertices);
  float3 vert = {0.0f,0.0f,0.0f};
  for (int i = 0; i < numVertices; ++i) {
    verticesTmp[i] = vert;
  }

  for (int i = 0; i < neighbors.size(); ++i) {
    std::set<int>::iterator it;
    for (it = neighbors[i].begin(); it != neighbors[i].end(); ++it) {
      double p[3];
      vertices->GetPoint(*it, p);
      verticesTmp[i] += make_float3(p[0],p[1],p[2]);
    }
  }

  for (int i = 0; i < verticesTmp.size(); ++i) {
    if (neighbors.size() > 5) {
      double p[3] = {verticesTmp[i].x/neighbors[i].size(),
		     verticesTmp[i].y/neighbors[i].size(),
		     verticesTmp[i].z/neighbors[i].size()};
      vertices->SetPoint(i,p);
    }
  }
}

bool cellOnInterface(vtkDataArray *data, int res[3], int i, int j, int k)
{
  int idx_left =   i-1 + j*res[0] +    k*res[0]*res[1];
  int idx_right =  i+1 + j*res[0] +    k*res[0]*res[1];
  int idx_bottom = i +  (j-1)*res[0] + k*res[0]*res[1];
  int idx_top =    i +  (j+1)*res[0] + k*res[0]*res[1];
  int idx_back =   i +   j*res[0] +   (k-1)*res[0]*res[1];
  int idx_front =  i +   j*res[0] +   (k+1)*res[0]*res[1];
  int idx = i + j*res[0] + k*res[0]*res[1];
  float f = data->GetComponent(idx,0);
  if (f > g_emf0 && f < g_emf1) {
    return true;
  }
  else if (f >= g_emf1) {
    if ((i-1 >= 0 && data->GetComponent(idx-1, 0) <= g_emf0) || 
	(i+1 < res[0] && data->GetComponent(idx+1, 0) <= g_emf0) || 
	(j-1 >= 0 && data->GetComponent(idx-res[0], 0) <= g_emf0) || 
	(j+1 < res[1] && data->GetComponent(idx+res[0], 0) <= g_emf0) || 
	(k-1 >= 0 && data->GetComponent(idx-res[0]*res[1], 0) <= g_emf0) || 
	(k+1 < res[2] && data->GetComponent(idx+res[0]*res[1], 0) <= g_emf0)) {
      return true;
    }
  }

  return false;
}

void generateSeedPoints(vtkRectilinearGrid *input,
			int refinement,
			vtkPoints *points,
			vtkIntArray *connectivity,
			vtkShortArray *coords)
{
  vtkDataArray *data =
    input->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  int inputRes[3];
  input->GetDimensions(inputRes);
  int cellRes[3] = {inputRes[0]-1, inputRes[1]-1, inputRes[2]-1};

  vtkDataArray *coordCenters[3];
  vtkDataArray *coordNodes[3];

  coordNodes[0] = input->GetXCoordinates();
  coordNodes[1] = input->GetYCoordinates();
  coordNodes[2] = input->GetZCoordinates();

  for (int c = 0; c < 3; ++c) {
    coordCenters[c] = vtkFloatArray::New();
    coordCenters[c]->SetNumberOfComponents(1);
    coordCenters[c]->SetNumberOfTuples(coordNodes[c]->GetNumberOfTuples()-1);
    for (int i = 0; i < coordCenters[c]->GetNumberOfTuples(); ++i) {
      coordCenters[c]->
	SetComponent(i,0,(coordNodes[c]->GetComponent(0,i) + 
			  coordNodes[c]->GetComponent(0,i+1))/2.0f);
    }
  }

  double bounds[6];
  input->GetBounds(bounds);
  int extent[6];
  input->GetExtent(extent);

  std::map<int3, int, bool(*)(const int3 &a, const int3 &b)> seedPos(compare_int3);
  seedPos.clear();
  int seedIdx = 0;

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

	  float gradf[3];
	  computeGradient(data, cellRes, i, j, k, coordCenters, gradf);

	  placeSeeds(points, cellCenter, cellSize, refinement, f, gradf,
		     bounds, i+extent[0], j+extent[2], k+extent[4], seedPos,
		     seedIdx);
	}
  	++idx;
      }
    }
  }

  connectivity->SetName("Connectivity");
  connectivity->SetNumberOfComponents(3);
  connectivity->SetNumberOfTuples(seedPos.size());

  coords->SetName("Coords");
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(seedPos.size());

  std::map<int3, int>::iterator it;
  for (it = seedPos.begin(); it != seedPos.end(); ++it) {

    int conn[3] = {-1,-1,-1};

    int x = it->first.x;
    int y = it->first.y;
    int z = it->first.z;

    int3 pos_xm = {x-1,y,z};
    if (seedPos.find(pos_xm) != seedPos.end()) {
      conn[0] = seedPos[pos_xm];
    }
    int3 pos_ym = {x,y-1,z};
    if (seedPos.find(pos_ym) != seedPos.end()) {
      conn[1] = seedPos[pos_ym];
    }
    int3 pos_zm = {x,y,z-1};
    if (seedPos.find(pos_zm) != seedPos.end()) {
      conn[2] = seedPos[pos_zm];
    }
    
    int seedIdx = it->second;
    connectivity->SetTuple3(seedIdx, conn[0], conn[1], conn[2]);
    coords->SetTuple3(seedIdx, x, y, z);
  }
}

float dot(float* a, float* b)
{ 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

float length(float* v)
{
  return sqrt(dot(v, v));
}

void normalize(float* v, float* n)
{
  float len = length(v);
  float invLen;

  if (len != 0.0f)
    invLen = 1.0f/len;
  else
    invLen = 0.0f;

  n[0] = v[0]*invLen;
  n[1] = v[1]*invLen;
  n[2] = v[2]*invLen;
}

void cross(float* a, float* b, float* c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

#define PI 3.14159265

void computeNormals(int nodeRes[3],
		    std::vector<float> &dx,
		    std::vector<float> &dy,
		    std::vector<float> &dz, 
		    vtkDataArray *f,
		    std::vector<float> &normals)
{
  // const double contact = 90;
  int i, j, k;
  float dfm1, dfm2;

  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

  for (k = 0; k < nodeRes[2]; k++) {
    int km = k - 1;
    int kp = k;
    if (km < 0) 
      km = 0;
    if (kp > cellRes[2]-1) 
      kp = cellRes[2]-1;

    float dzc = (dz[km] + dz[kp])*0.5f;
    
    for (j = 0; j < nodeRes[1]; j++) {
      int jm = j - 1;
      int jp = j;
      if (jm < 0) 
	jm = 0;
      if (jp > cellRes[1]-1) 
      	jp = cellRes[1]-1;

      float dyc = (dy[jm] + dy[jp])*0.5f;

      for (i = 0; i < nodeRes[0]; i++) {
	int im = i - 1;
	int ip = i;
	if (im < 0) 
	  im = 0;
	if (ip > cellRes[0]-1) 
	  ip = cellRes[0]-1;
	
	float dxc = (dx[im] + dx[ip])*0.5f;

	float fs[8] = {f->GetComponent(im + jm*cellRes[0] + km*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(ip + jm*cellRes[0] + km*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(im + jp*cellRes[0] + km*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(ip + jp*cellRes[0] + km*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(im + jm*cellRes[0] + kp*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(ip + jm*cellRes[0] + kp*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(im + jp*cellRes[0] + kp*cellRes[0]*cellRes[1], 0),
		       f->GetComponent(ip + jp*cellRes[0] + kp*cellRes[0]*cellRes[1], 0)};


	dfm1 = (fs[7] - fs[6])*dz[km] + (fs[3] - fs[2])*dz[kp];
	dfm2 = (fs[5] - fs[4])*dz[km] + (fs[1] - fs[0])*dz[kp];
	float nx = 0.25f*(dfm1*dy[j]+dfm2*dy[jp]) / (dxc*dyc*dzc);

	dfm1 = (fs[7] - fs[5])*dz[km] + (fs[3] - fs[1])*dz[kp];
	dfm2 = (fs[6] - fs[4])*dz[km] + (fs[2] - fs[0])*dz[kp];
	float ny = 0.25f*(dfm1*dx[i]+dfm2*dx[ip]) / (dxc*dyc*dzc);

	dfm1 = (fs[7] - fs[3])*dy[jm] + (fs[5] - fs[1])*dy[jp];
	dfm2 = (fs[6] - fs[2])*dy[jm] + (fs[4] - fs[0])*dy[jp];
	float nz = 0.25f*(dfm1*dx[i]+dfm2*dx[ip]) / (dxc*dyc*dzc);

	int offset = i+j*nodeRes[0]+k*nodeRes[0]*nodeRes[1];

	normals[offset*3+0] = -nx;// normal points from f outwards
	normals[offset*3+1] = -ny;
	normals[offset*3+2] = -nz;
      }
    }
  }
}

float computeLstar(float f, float n[3], float d[3])
{
  const int MAXITER = 100;
  int ih, i1, i2, i3;
  float d1, d2, d3;
  float n1, n2, n3;
  float nd[3], nd1, nd2, nd3, ndsum;
  float volume;
  float li, lii, liii, liv, lv;
  float la, lb, ll, sla, slb, sumla, sumlb, dlstar;
  int niter;
  float d2d3, n2rd3, n2n3;
  float ti, tii, tiii, tiv, tv;
  float lstar;
  float epsf = 0.001f; // error (?)

  niter = 0;
  nd[0] = fabs(n[0]*d[0]);
  nd[1] = fabs(n[1]*d[1]);
  nd[2] = fabs(n[2]*d[2]);
  i1 = 0; // indices decremented by 1
  i2 = 1;
  i3 = 2;

  // [i3] < [i2] < [i1] (?)
  if (nd[0] < nd[1]) {
    i1 = 1;
    i2 = 0;
  }
  if (nd[i2] < nd[2]) {
    i3 = i2;
    i2 = 2;
  }
  if (nd[i1] < nd[i2]) {
    ih = i1;
    i1 = i2;
    i2 = ih;
  }

  d1 = d[i1];
  d2 = d[i2];
  d3 = d[i3];

  n1 = fabs(n[i1]);
  n2 = fabs(n[i2]);
  n3 = fabs(n[i3]);

  nd1 = nd[i1];
  nd2 = nd[i2];
  nd3 = nd[i3];

  ndsum = nd1 + nd2 + nd3;
  volume = d1 * d2 * d3; // so d is box dimension ?

  d2d3 = d2 * d3;
  n2rd3 = n2 / d3;
  n2n3 = n2 * n3;

  if (f < g_emf0)
    return 0.0f;
  if (f > g_emf1)
    return ndsum;

  li = nd1;
  lii = nd1 + nd3;
  liii = nd1 + nd2;
  liv = liii + nd3;
  lv = liv + nd1;

  dlstar = 0.0f;
  lstar = 0.5f * liv;


  sumla = 0.0f;
  sumlb = 0.5f * volume * n1;

  while (1) {
    if (fabs((sumlb-sumla)/volume/n1 - f) < epsf || niter > MAXITER)
      break;
    niter = niter + 1;
    la = lstar;
    lb = lstar + nd1;

    //-------------------------------------------------------------------------
    // calculation of ll, sla, slb
    // Bereich 1
    if (la >= 0.0f && la <= li)
      {
	sla = 0.0f;
	sumla = 0.0f;
      }
    else if (lb >= 0.0f && lb <= li)
      {
	slb = 0.0f;
	sumlb = 0.0f;
      }
    // Bereich 2
    if (la >= li && la <= lii)
      {
	ll = la - li;
	sla = ll*ll / (2.0f*n2n3);
	sumla = ll*ll*ll/ (6.0f*n2n3);
      }
    else if (lb >= li && lb <= lii)
      {
	ll = lb - li;
	slb = ll*ll / (2.0f*n2n3);
	sumlb = ll*ll*ll/ (6.0f*n2n3);
      }
    // Bereich 3
    if (la >= lii && la <= liii)
      {
	ll = la - lii;
	sla = nd3 / n2rd3 / 2.0f + ll / n2rd3;
	sumla = (3.0f*ll*(nd3 + ll) + nd3*nd3) / n2rd3 / 6.0f;
      }
    else if (lb >= lii && lb <= liii)
      {
	ll = lb - lii;
	slb = nd3 / n2rd3 / 2.0f + ll / n2rd3;
	sumlb = (3.0f*ll*(nd3 + ll) + nd3*nd3) / n2rd3 / 6.0f;
      }
    // Bereich 4
    if (la >= liii && la <= liv)
      {
	ll = liv - la;
	sla = d2d3 - 0.5f*ll*ll / n2n3;
	sumla = (3.0f*nd2*nd3*nd3 + 3.0f*nd2*nd2*nd3 -
		 6.0f*nd2*nd3*ll + ll*ll*ll) / n2n3 / 6.0f;
      }
    else if (lb >= liii && lb <= liv)
      {
	ll = liv - lb;
	slb = d2d3 - 0.5f*ll*ll / n2n3;
	sumlb = (3.0f*nd2*nd3*nd3 + 3.0f*nd2*nd2*nd3 -
		 6.0f*nd2*nd3*ll + ll*ll*ll) / n2n3 / 6.0f;
      }
    // Bereich 5
    if (la >= liv && la <= lv)
      {
	ll = la - liv;
	sla = d2d3;
	sumla = d2d3*(ll + 0.5f*nd2 + 0.5f*nd3);
      }
    else if (lb >= liv && lb <= lv)
      {
	ll = lb - liv;
	slb = d2d3;
	sumlb = d2d3*(ll + 0.5f*nd2 + 0.5f*nd3);
      }
    dlstar = (sumlb - sumla - f*volume*n1) / (slb - sla);
    lstar = lstar - dlstar;
    lstar = std::max(lstar, 0.0f);
    lstar = std::min(lstar, liv);
  }

  return lstar;
}

void computeL(int cellRes[3],
	      std::vector<float> &dx,
	      std::vector<float> &dy,
	      std::vector<float> &dz, 
	      vtkDataArray *f,
	      std::vector<float> &normals,
	      std::vector<float> &lstar,
	      std::vector<float> &normalsInt)
{
  //  int no = 0;
  int w, h, d;
  w = cellRes[0];
  h = cellRes[1];
  d = cellRes[2];

  for (int k = 0; k < d; k++) {
    for (int j = 0; j < h; j++) {
      for (int i = 0; i < w; i++) {
	int fo = i + j*w + k*w*h;

	// The correct normals vector is computed as an average of 
	// 8 corners;
	int n0 = i   + (j  )*(w+1) + (k  )*(w+1)*(h+1);
	int n1 = i+1 + (j  )*(w+1) + (k  )*(w+1)*(h+1);
	int n2 = i   + (j+1)*(w+1) + (k  )*(w+1)*(h+1);
	int n3 = i+1 + (j+1)*(w+1) + (k  )*(w+1)*(h+1);
	int n4 = i   + (j  )*(w+1) + (k+1)*(w+1)*(h+1);
	int n5 = i+1 + (j  )*(w+1) + (k+1)*(w+1)*(h+1);
	int n6 = i   + (j+1)*(w+1) + (k+1)*(w+1)*(h+1);
	int n7 = i+1 + (j+1)*(w+1) + (k+1)*(w+1)*(h+1);

	float ns[8][3] = {{normals[n0*3+0], normals[n0*3+1], normals[n0*3+2]},
			  {normals[n1*3+0], normals[n1*3+1], normals[n1*3+2]},
			  {normals[n2*3+0], normals[n2*3+1], normals[n2*3+2]},
			  {normals[n3*3+0], normals[n3*3+1], normals[n3*3+2]},
			  {normals[n4*3+0], normals[n4*3+1], normals[n4*3+2]},
			  {normals[n5*3+0], normals[n5*3+1], normals[n5*3+2]},
			  {normals[n6*3+0], normals[n6*3+1], normals[n6*3+2]},
			  {normals[n7*3+0], normals[n7*3+1], normals[n7*3+2]}};

	float n[3] = {0.0f, 0.0f, 0.0f};

	for (int l = 0; l < 8; l++)
	  {
	    n[0] += ns[l][0];
	    n[1] += ns[l][1];
	    n[2] += ns[l][2];
	  } 
	float len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	if (len)
	  {
	    n[0] /= len;
	    n[1] /= len;
	    n[2] /= len;
	  }

	normalsInt[fo*3+0] = n[0];
	normalsInt[fo*3+1] = n[1];
	normalsInt[fo*3+2] = n[2];

	float dd[3] = {dx[i], dy[j], dz[k]};
	
	if (f->GetComponent(fo,0) > g_emf0 && f->GetComponent(fo,0) < g_emf1) {
	  lstar[fo] = computeLstar(f->GetComponent(fo,0), n, dd);
	}
	else if (f->GetComponent(fo,0) >= g_emf1) {
	  if (f->GetComponent(fo-1,0) < g_emf0 || f->GetComponent(fo+1,0) < g_emf0 || 
	      f->GetComponent(fo-w,0) < g_emf0 || f->GetComponent(fo+w,0) < g_emf0 || 
	      f->GetComponent(fo-w*h,0) < g_emf0 || f->GetComponent(fo+w*h,0) < g_emf0) {
	    lstar[fo] = computeLstar(f->GetComponent(fo,0), n, dd);
	  }
	}
	else { 
	  lstar[fo] = 0.0f;
	}
      }
    }
  }
}

void generateSeedPointsPLIC(vtkRectilinearGrid *vofGrid,
			    int refinement,
			    vtkPoints *points,
			    vtkIntArray *connectivity,
			    vtkShortArray *coords,
			    int globalExtent[6],
			    int numGhostLevels)
{  
  vtkDataArray *vofArray =
    vofGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  vtkDataArray *coordNodes[3] = {vofGrid->GetXCoordinates(),
				 vofGrid->GetYCoordinates(),
				 vofGrid->GetZCoordinates()};
  int nodeRes[3] = {coordNodes[0]->GetNumberOfTuples(),
		    coordNodes[1]->GetNumberOfTuples(),
		    coordNodes[2]->GetNumberOfTuples()};
  int cellRes[3] = {nodeRes[0]-1,nodeRes[1]-1,nodeRes[2]-1};
  
  std::vector<std::vector<float> > dx(3);
  dx[0].resize(cellRes[0]);
  dx[1].resize(cellRes[1]);
  dx[2].resize(cellRes[2]);

  for (int c = 0; c < 3; ++c) {
    for (int i = 0; i < cellRes[c]; ++i) {
      dx[c][i] = coordNodes[c]->GetComponent(i+1,0) - coordNodes[c]->GetComponent(i,0);
    }
  }

  std::vector<float> normals;  
  normals.resize(nodeRes[0]*nodeRes[1]*nodeRes[2]*3);

  computeNormals(nodeRes, dx[0], dx[1], dx[2], vofArray, normals);

  std::vector<float> lstar(cellRes[0]*cellRes[1]*cellRes[2]);
  std::vector<float> normalsInt(cellRes[0]*cellRes[1]*cellRes[2]*3);
  computeL(cellRes, dx[0], dx[1], dx[2], vofArray, normals, lstar, normalsInt);

  //--------------
  vtkDataArray *data =
    vofGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  int inputRes[3];
  vofGrid->GetDimensions(inputRes);
  vtkDataArray *coordCenters[3];

  for (int c = 0; c < 3; ++c) {
    coordCenters[c] = vtkFloatArray::New();
    coordCenters[c]->SetNumberOfComponents(1);
    coordCenters[c]->SetNumberOfTuples(coordNodes[c]->GetNumberOfTuples()-1);
    for (int i = 0; i < coordCenters[c]->GetNumberOfTuples(); ++i) {
      coordCenters[c]->
	SetComponent(i,0,(coordNodes[c]->GetComponent(0,i) + 
			  coordNodes[c]->GetComponent(0,i+1))/2.0f);
    }
  }

  double bounds[6];
  vofGrid->GetBounds(bounds);
  int extent[6];
  vofGrid->GetExtent(extent);

  std::map<int3, int, bool(*)(const int3 &a, const int3 &b)> seedPos(compare_int3);
  seedPos.clear();
  int seedIdx = 0;

  //---------------------------------------------------------------------------
  // populate the grid with seed points
  // int idx = 0;
  float cellCenter[3];
  float cellSize[3];

  int imin = extent[0] > globalExtent[0] ? numGhostLevels : 0;
  int imax = extent[1] < globalExtent[1] ? cellRes[0]-numGhostLevels : cellRes[0];
  int jmin = extent[2] > globalExtent[2] ? numGhostLevels : 0;
  int jmax = extent[3] < globalExtent[3] ? cellRes[1]-numGhostLevels : cellRes[1];
  int kmin = extent[4] > globalExtent[4] ? numGhostLevels : 0;
  int kmax = extent[5] < globalExtent[5] ? cellRes[2]-numGhostLevels : cellRes[2];

  int kcur = kmin;
  for (int k = kmin; k < kmax; ++k) {
    cellCenter[2] = coordCenters[2]->GetComponent(kcur,0);
    cellSize[2] = coordNodes[2]->GetComponent(kcur+1,0) - coordNodes[2]->GetComponent(kcur,0);

    int jcur = jmin;
    for (int j = jmin; j < jmax; ++j) {      
      cellCenter[1] = coordCenters[1]->GetComponent(jcur,0);
      cellSize[1] = coordNodes[1]->GetComponent(jcur+1,0) - coordNodes[1]->GetComponent(jcur,0);

      int icur = imin;
      for (int i = imin; i < imax; ++i) {	
	cellCenter[0] = coordCenters[0]->GetComponent(icur,0);
	cellSize[0] = coordNodes[0]->GetComponent(icur+1,0) - coordNodes[0]->GetComponent(icur,0);

	int idx = i + j*cellRes[0] + k*cellRes[0]*cellRes[1];
	float f = data->GetComponent(0,idx);
	if (f > g_emf0) {

	  placeSeedsPLIC(points, cellCenter, cellSize, refinement, cellRes, f, lstar, normalsInt,
			 bounds, i, j, k, idx, seedPos, seedIdx);
	}
	++icur;
      }
      ++jcur;
    }
    ++kcur;
  }

  connectivity->SetName("Connectivity");
  connectivity->SetNumberOfComponents(3);
  connectivity->SetNumberOfTuples(seedPos.size());

  coords->SetName("Coords");
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(seedPos.size());
  
  std::map<int3, int>::iterator it;
  for (it = seedPos.begin(); it != seedPos.end(); ++it) {

    int conn[3] = {-1,-1,-1};

    int x = it->first.x;
    int y = it->first.y;
    int z = it->first.z;

    int3 pos_xm = {x-1,y,z};
    if (seedPos.find(pos_xm) != seedPos.end()) {
      conn[0] = seedPos[pos_xm];
    }
    int3 pos_ym = {x,y-1,z};
    if (seedPos.find(pos_ym) != seedPos.end()) {
      conn[1] = seedPos[pos_ym];
    }
    int3 pos_zm = {x,y,z-1};
    if (seedPos.find(pos_zm) != seedPos.end()) {
      conn[2] = seedPos[pos_zm];
    }
    
    int seedIdx = it->second;
    connectivity->SetTuple3(seedIdx, conn[0], conn[1], conn[2]);
    coords->SetTuple3(seedIdx, x, y, z);
  }
}


// void advectParticles(vtkRectilinearGrid *vofGrid[2],
// 		     vtkRectilinearGrid *velocityGrid[2],
// 		     std::vector<float4> &particles,
// 		     const float deltaT)
// {
//   int nodeRes[3];
//   vofGrid[0]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
//   vtkDataArray *velocityArray0 = velocityGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *velocityArray1 = velocityGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *vofArray0 = vofGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
//   vtkDataArray *vofArray1 = vofGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

//   int numSubsteps = 10;
  
//   for (int i = 0; i < particles.size(); ++i) {

//     float4 p = particles[i];
//     float deltaTSub = deltaT/numSubsteps;

//     for (int j = 0; j < numSubsteps; ++j) {

//       double x[3] = {p.x, p.y, p.z};
//       int ijk[3];
//       double pcoords[3];
      
//       velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
//       float4 velocity0 = make_float4(interpolateVec(velocityArray0, cellRes, ijk, pcoords),0.0f);
//       float4 velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);

//       float a = float(j)/numSubsteps;

//       float4 velocity = (1.0f-a)*velocity0 + a*velocity1;

//       if (j == numSubsteps-1) {
// 	deltaTSub = deltaT - j*deltaTSub;
//       }
//       p += velocity * deltaTSub;
//     }

//     particles[i] = p;

//     double x1[3] = {particles[i].x, particles[i].y, particles[i].z};
//     int ijk[3];
//     double pcoords[3];
//     vofGrid[1]->ComputeStructuredCoordinates(x1, ijk, pcoords);
//     int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
//     float f = vofArray1->GetComponent(idx, 0);
//     if (f <= g_emf0) {
//       particles[i].w = 0.0f;
//     }   
//   }
// }

// 2615
// Heun's
// void advectParticles(vtkRectilinearGrid *vofGrid[2],
// 		     vtkRectilinearGrid *velocityGrid[2],
// 		     std::vector<float4> &particles,
// 		     const float deltaT)
// {
//   int nodeRes[3];
//   vofGrid[0]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
//   vtkDataArray *velocityArray0 = velocityGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *velocityArray1 = velocityGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *vofArray0 = vofGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
//   vtkDataArray *vofArray1 = vofGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
  
//   for (int i = 0; i < particles.size(); ++i) {

//     double x[3] = {particles[i].x, particles[i].y, particles[i].z};
//     int ijk[3];
//     double pcoords[3];
//     vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
//     float4 velocity0 = make_float4(interpolateVec(velocityArray0, cellRes, ijk, pcoords),0.0f);

//     float4 p1 = particles[i] + velocity0*deltaT;
//     x[0] = p1.x;
//     x[1] = p1.y;
//     x[2] = p1.z;    
//     vofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
//     float4 velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
//     particles[i] += 0.5f*(velocity0 + velocity1)*deltaT;

//     // particles[i] += velocity0*deltaT;

//     x[0] = particles[i].x;
//     x[1] = particles[i].y;
//     x[2] = particles[i].z;
//     vofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
//     int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
//     float f = vofArray1->GetComponent(idx, 0);
//     if (f <= g_emf0) {
//       particles[i].w = 0.0f;
//     }   
//   }
// }

// //2100
// void advectParticles(vtkRectilinearGrid *vofGrid[2],
// 		     vtkRectilinearGrid *velocityGrid[2],
// 		     std::vector<float4> &particles,
// 		     const float deltaT)
// {
//   int nodeRes[3];
//   vofGrid[0]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
//   vtkDataArray *velocityArray = velocityGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *vofArray = vofGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

//   // ---------------------------
//   vtkDataArray *coordCenters[3];
//   vtkDataArray *coordNodes[3];

//   coordNodes[0] = vofGrid[0]->GetXCoordinates();
//   coordNodes[1] = vofGrid[0]->GetYCoordinates();
//   coordNodes[2] = vofGrid[0]->GetZCoordinates();

//   for (int c = 0; c < 3; ++c) {
//     coordCenters[c] = vtkFloatArray::New();
//     coordCenters[c]->SetNumberOfComponents(1);
//     coordCenters[c]->SetNumberOfTuples(coordNodes[c]->GetNumberOfTuples()-1);
//     for (int i = 0; i < coordCenters[c]->GetNumberOfTuples(); ++i) {
//       coordCenters[c]->
//   	SetComponent(i,0,(coordNodes[c]->GetComponent(0,i) + 
//   			  coordNodes[c]->GetComponent(0,i+1))/2.0f);
//     }
//   }
//   float dx = coordCenters[0]->GetComponent(1,0) - coordCenters[0]->GetComponent(0,0);
//   float dy = coordCenters[1]->GetComponent(1,0) - coordCenters[1]->GetComponent(0,0);
//   float dz = coordCenters[2]->GetComponent(1,0) - coordCenters[2]->GetComponent(0,0);  
//   // ---------------------------

//   std::vector<float4>::iterator it;
//   for (it = particles.begin(); it != particles.end(); ++it) {

//     if (it->w != 0.0f) {
//       double x[3] = {it->x, it->y, it->z};
//       int ijk[3];
//       double pcoords[3];
//       int particleInsideGrid = vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
//       int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
//       float f = vofArray->GetComponent(idx, 0);
      
//       if (particleInsideGrid) {
// 	// ---------------------------
// 	int iter = 0;
// 	do {
// 	  iter++;
// 	  float grad[3];
// 	  computeGradient(vofArray,  cellRes, ijk[0], ijk[1], ijk[2], coordCenters, grad);
// 	  // computeGradient(vofGrid[0], vofArray, cellRes, ijk, coordCenters, pcoords, grad);
// 	  float3 gr = make_float3(grad[0],grad[1],grad[2]);
// 	  if (length(gr) > 0.0f)
// 	    gr = normalize(gr);
// 	  x[0] += gr.x*dx;
// 	  x[1] += gr.y*dy;
// 	  x[2] += gr.z*dz;
	  
// 	  particleInsideGrid = vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
	  
// 	  idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
// 	  f = vofArray->GetComponent(idx, 0);	  
// 	}  while (f <= g_emf0 && iter < 10);
// 	if (f > g_emf0) {
// 	  *it = make_float4(x[0],x[1],x[2],1.0f);
// 	}
// 	else {
// 	  it->w = 0.0f;
// 	  continue;
// 	}
// 	// ---------------------------
// 	float4 velocity = make_float4(interpolateVec(velocityArray, cellRes, ijk, pcoords),0.0f);
// 	*it = *it + velocity*deltaT;
//       }
//     }
//   }
// }


// // dirty hack - sending particles to the nearest cell with f > 0
// //2100
// void advectParticles(vtkRectilinearGrid *vofGrid[2],
// 		     vtkRectilinearGrid *velocityGrid[2],
// 		     std::vector<float4> &particles,
// 		     const float deltaT)
// {
//   int nodeRes[3];
//   vofGrid[0]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
//   vtkDataArray *velocityArray = velocityGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *vofArray1 = vofGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

//   // ---------------------------
//   vtkDataArray *coordCenters[3];
//   vtkDataArray *coordNodes[3];

//   coordNodes[0] = vofGrid[0]->GetXCoordinates();
//   coordNodes[1] = vofGrid[0]->GetYCoordinates();
//   coordNodes[2] = vofGrid[0]->GetZCoordinates();

//   for (int c = 0; c < 3; ++c) {
//     coordCenters[c] = vtkFloatArray::New();
//     coordCenters[c]->SetNumberOfComponents(1);
//     coordCenters[c]->SetNumberOfTuples(coordNodes[c]->GetNumberOfTuples()-1);
//     for (int i = 0; i < coordCenters[c]->GetNumberOfTuples(); ++i) {
//       coordCenters[c]->
//   	SetComponent(i,0,(coordNodes[c]->GetComponent(0,i) + 
//   			  coordNodes[c]->GetComponent(0,i+1))/2.0f);
//     }
//   }
//   float dx = coordCenters[0]->GetComponent(1,0) - coordCenters[0]->GetComponent(0,0);
//   float dy = coordCenters[1]->GetComponent(1,0) - coordCenters[1]->GetComponent(0,0);
//   float dz = coordCenters[2]->GetComponent(1,0) - coordCenters[2]->GetComponent(0,0);  
//   // ---------------------------

//   std::vector<float4>::iterator it;

//   int numStranded = 0;
//   for (it = particles.begin(); it != particles.end(); ++it) {
//     if (it->w == 0.0f)
//       numStranded++;
//   }
//   std::cout << "numStranded before = " << numStranded << std::endl;
  
//   for (it = particles.begin(); it != particles.end(); ++it) {

//     if (it->w != 0.0f) {
      
//       double x[3] = {it->x, it->y, it->z};
//       int ijk[3];
//       double pcoords[3];
//       vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
      
//       float4 velocity = make_float4(interpolateVec(velocityArray, cellRes, ijk, pcoords),0.0f);
//       *it = *it + velocity*deltaT;
//       x[0] = it->x;
//       x[1] = it->y;
//       x[2] = it->z;
//       int particleInsideGrid = vofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
//       int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
//       float f = vofArray1->GetComponent(idx, 0);
     
//       if (particleInsideGrid) {

// 	if (f <= g_emf0) {
// 	  // find nearest cell with f > g_emf0
// 	  int tol = 25;
// 	  float minDist = std::numeric_limits<float>::max();
// 	  int ci, cj, ck;
// 	  for (int rk = -tol; rk <= tol; ++rk) {
// 	    if (rk+ijk[2] < 0 || rk+ijk[2] >= cellRes[2])
// 	      continue;
// 	    for (int rj = -tol; rj <= tol; ++rj) {
// 	      if (rj+ijk[1] < 0 || rj+ijk[1] >= cellRes[1])
// 		continue;
// 	      for (int ri = -tol; ri <= tol; ++ri) {
// 		if (ri+ijk[0] < 0 || ri+ijk[0] >= cellRes[0])
// 		  continue;

// 		int idx = ri+ijk[0] + (rj+ijk[1])*cellRes[0] + (rk+ijk[2])*cellRes[0]*cellRes[1];
// 		float f = vofArray1->GetComponent(idx, 0);
// 		if (f > g_emf0) {
// 		  float dist = ri*ri + rj*rj + rk*rk;
// 		  if (dist < minDist) {
// 		    minDist = dist;
// 		    ci = ri;
// 		    cj = rj;		
// 		    ck = rk;
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	  if (minDist < std::numeric_limits<float>::max()) {
// 	    *it = make_float4(x[0] + ci*dx,
// 			      x[1] + cj*dy,
// 			      x[2] + ck*dz,
// 			      1.0f);
// 	  }
// 	  else {
// 	    it->w = 0.0f;
// 	  }
// 	}
// 	// int iter = 0;
// 	// while (iter < 10 && f <= g_emf0) {
// 	//   iter++;
// 	//   float grad[3];
// 	//   computeGradient(vofArray1,  cellRes, ijk[0], ijk[1], ijk[2], coordCenters, grad);
// 	//   float3 gr = make_float3(grad[0],grad[1],grad[2]);
// 	//   if (length(gr) > 0.0f)
// 	//     gr = normalize(gr);
// 	//   x[0] += gr.x*dx;
// 	//   x[1] += gr.y*dy;
// 	//   x[2] += gr.z*dz;

// 	//   vofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);	  
// 	//   idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
// 	//   f = vofArray1->GetComponent(idx, 0);	  
// 	// }
// 	// if (f > g_emf0) {
// 	//   *it = make_float4(x[0],x[1],x[2],1.0f);
// 	// }
// 	// else {
// 	//   it->w = 0.0f;
// 	// }
// 	// ---------------------------
//       }      
//       else {
// 	it->w = 0.0f;
//       }
//     }
//   }
//   numStranded = 0;
//   for (it = particles.begin(); it != particles.end(); ++it) {
//     if (it->w == 0.0f)
//       numStranded++;
//   }
//   std::cout << "numStranded after = " << numStranded << std::endl;

// }

// iterative
// https://en.wikipedia.org/wiki/Trapezoidal_rule_%28differential_equations%29
// 2100
void advectParticles(vtkRectilinearGrid *vofGrid[2],
		     vtkRectilinearGrid *velocityGrid[2],
		     std::vector<float4> &particles,
		     const float deltaT)
{
  int nodeRes[3];
  vofGrid[0]->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  vtkDataArray *velocityArray0 = velocityGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
  vtkDataArray *velocityArray1 = velocityGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
  vtkDataArray *vofArray1 = vofGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

  std::vector<float4>::iterator it;
  for (it = particles.begin(); it != particles.end(); ++it) {

    if (it->w != 0.0f) {
      
      double x[3];
      int ijk[3];
      double pcoords[3];
      const int maxNumIter = 20;
      
      float4 pos0 = *it;
      x[0] = pos0.x;
      x[1] = pos0.y;
      x[2] = pos0.z;        
      vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
      float4 velocity0 = make_float4(interpolateVec(velocityArray0, cellRes, ijk, pcoords),0.0f);
      float4 velocity = velocity0;
      
      for (int i = 0; i < maxNumIter; ++i) {

	float4 pos1 = pos0 + velocity * deltaT;
	x[0] = pos1.x;
	x[1] = pos1.y;
	x[2] = pos1.z;
	velocityGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
	float4 velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
	
	velocity = (velocity0 + velocity1)/2.0f;
      }
      
      *it = *it + velocity*deltaT;

      x[0] = it->x;
      x[1] = it->y;
      x[2] = it->z;
      int particleInsideGrid = vofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
      int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      float f = vofArray1->GetComponent(idx, 0);
            
      if (particleInsideGrid) {
	if (f <= g_emf0) {
	  it->w = 0.0f;
	}
      }
    }
  }
}

// // basic
// //2100
// void advectParticles(vtkRectilinearGrid *vofGrid[2],
// 		     vtkRectilinearGrid *velocityGrid[2],
// 		     std::vector<float4> &particles,
// 		     const float deltaT)
// {
//   int nodeRes[3];
//   vofGrid[0]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
//   vtkDataArray *velocityArray = velocityGrid[0]->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
//   vtkDataArray *vofArray1 = vofGrid[1]->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

//   std::vector<float4>::iterator it;
//   for (it = particles.begin(); it != particles.end(); ++it) {

//     if (it->w != 0.0f) {
      
//       double x[3] = {it->x, it->y, it->z};
//       int ijk[3];
//       double pcoords[3];
//       vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
//       float4 velocity = make_float4(interpolateVec(velocityArray, cellRes, ijk, pcoords),0.0f);
//       *it = *it + velocity*deltaT;

//       x[0] = it->x;
//       x[1] = it->y;
//       x[2] = it->z;
//       int particleInsideGrid = vofGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
//       int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
//       float f = vofArray1->GetComponent(idx, 0);
            
//       if (particleInsideGrid) {
// 	if (f <= g_emf0) {
// 	  it->w = 0.0f;
// 	}
//       }
//       else {
// 	it->w = 0.0f;
//       }
//     }
//   }
// }

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

int outOfBounds(const float4 particle, const double bounds[6], const double globalBounds[6])
{
  if (particle.x < bounds[0] && bounds[0] != globalBounds[0]) return 0;
  if (particle.x > bounds[1] && bounds[1] != globalBounds[1]) return 1;
  if (particle.y < bounds[2] && bounds[2] != globalBounds[2]) return 2;
  if (particle.y > bounds[3] && bounds[3] != globalBounds[3]) return 3;
  if (particle.z < bounds[4] && bounds[4] != globalBounds[4]) return 4;
  if (particle.z > bounds[5] && bounds[5] != globalBounds[5]) return 5;

  return -1;
}

int withinBounds(const float4 particle, const double bounds[6])
{
  if (particle.x < bounds[0]) return 0;
  if (particle.x > bounds[1]) return 0;
  if (particle.y < bounds[2]) return 0;
  if (particle.y > bounds[3]) return 0;
  if (particle.z < bounds[4]) return 0;
  if (particle.z > bounds[5]) return 0;

  return 1;
}

void prepareLabelsToSend(std::vector<std::vector<int> > &NeighborProcesses,
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

void setupNodeOffsets(vtkPoints *points, vtkIntArray *connectivity,
		      float co[3][12], float coc[3][3],
		      float po[3][12], float poc[3][3])
{
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

  const float coTmp[3][12] = {{0,0,0, 0,0,1, 0,1,1, 0,1,0},
			      {0,0,0, 1,0,0, 1,0,1, 0,0,1},
			      {0,0,0, 0,1,0, 1,1,0, 1,0,0}};
  const float cocTmp[3][3] = {{0,0.5f,0.5f},
			      {0.5f,0,0.5f},
			      {0.5f,0.5f,0}};
  
  const float poTmp[3][12] = {{-cs2[0],-cs2[1],-cs2[2], 
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
			       cs2[0],-cs2[1],-cs2[2]}};
  const float pocTmp[3][3] = {{-cs2[0],0,0},
			      {0,-cs2[1],0},
			      {0,0,-cs2[2]}};

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 12; ++j) {
      co[i][j] = coTmp[i][j];
      po[i][j] = poTmp[i][j];
    }
    for (int j = 0; j < 3; ++j) {
      coc[i][j] = cocTmp[i][j];
      poc[i][j] = pocTmp[i][j];
    }
  }
}

void generateBoundaries(vtkPoints *points,
			vtkFloatArray *labels,
			vtkIntArray *connectivity,
			vtkShortArray *coords,
			vtkPolyData *boundaries)
{
  float co[3][12];
  float coc[3][3];  
  float po[3][12];
  float poc[3][3];

  setupNodeOffsets(points, connectivity, co, coc, po, poc);

  vtkFloatArray *labelsBack = vtkFloatArray::New();
  labelsBack->SetNumberOfComponents(1);
  labelsBack->SetName("BackLabels");
  vtkFloatArray *labelsFront = vtkFloatArray::New();
  labelsFront->SetNumberOfComponents(1);
  labelsFront->SetName("FrontLabels");
  
  std::vector<float3> ivertices;
  ivertices.clear();
  std::vector<float3> vertices;
  vertices.clear();

  int vidx = 0;

  const int numPoints = points->GetNumberOfPoints();
  for (int i = 0; i < numPoints; ++i) {

    int conn[3] = {connectivity->GetComponent(i, 0),
  		   connectivity->GetComponent(i, 1),
  		   connectivity->GetComponent(i, 2)};

    float l0 = labels->GetValue(i);
    float c0[3] = {coords->GetComponent(i, 0),
		   coords->GetComponent(i, 1),
		   coords->GetComponent(i, 2)};

    double p0[3];
    points->GetPoint(i, p0);

    for (int j = 0; j < 3; ++j) {
      if (conn[j] > -1) {

  	float l1 = labels->GetValue(conn[j]);
	
  	double p1[3];
  	points->GetPoint(conn[j], p1);
      
  	if (l0 != l1) {

  	  labelsBack->InsertNextValue(l0);
  	  labelsBack->InsertNextValue(l0);
  	  labelsBack->InsertNextValue(l0);
  	  labelsBack->InsertNextValue(l0);
  	  labelsFront->InsertNextValue(l1);
  	  labelsFront->InsertNextValue(l1);
  	  labelsFront->InsertNextValue(l1);
  	  labelsFront->InsertNextValue(l1);

  	  float3 verts[5];
  	  float3 iverts[5];

	  verts[4] = make_float3(p0[0]+poc[j][0], 
				 p0[1]+poc[j][1], 
				 p0[2]+poc[j][2]);
	  iverts[4] = make_float3(c0[0]+coc[j][0], 
				  c0[1]+coc[j][1], 
				  c0[2]+coc[j][2]);

  	  for (int k = 0; k < 4; ++k) {

  	    float3 vertex = {p0[0]+po[j][k*3+0], 
			     p0[1]+po[j][k*3+1], 
			     p0[2]+po[j][k*3+2]};
  	    verts[k] = vertex;

  	    float3 ivertex = {c0[0]+co[j][k*3+0], 
			      c0[1]+co[j][k*3+1], 
			      c0[2]+co[j][k*3+2]};
  	    iverts[k] = ivertex;
  	  }
  	  vertices.push_back(verts[0]);
  	  vertices.push_back(verts[1]);
	  vertices.push_back(verts[4]);
  	  vertices.push_back(verts[1]);
  	  vertices.push_back(verts[2]);
	  vertices.push_back(verts[4]);
  	  vertices.push_back(verts[2]);
  	  vertices.push_back(verts[3]);
	  vertices.push_back(verts[4]);
	  vertices.push_back(verts[3]);
	  vertices.push_back(verts[0]);
	  vertices.push_back(verts[4]);
	  
  	  ivertices.push_back(iverts[0]);
  	  ivertices.push_back(iverts[1]);
	  ivertices.push_back(iverts[4]);
  	  ivertices.push_back(iverts[1]);
  	  ivertices.push_back(iverts[2]);
	  ivertices.push_back(iverts[4]);
  	  ivertices.push_back(iverts[2]);
  	  ivertices.push_back(iverts[3]);
	  ivertices.push_back(iverts[4]);
  	  ivertices.push_back(iverts[3]);
  	  ivertices.push_back(iverts[0]);
	  ivertices.push_back(iverts[4]);

	  vidx += 12;
  	}
      }
    }
  }

  std::vector<int> indices;
  std::vector<float3> mergedVertices;

  mergeTriangles(vertices, ivertices, indices, mergedVertices);
  
  // for (int i = 0; i < 10; ++i)
  // smoothSurface(mergedVertices, indices);
  
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

  boundaries->SetPoints(outputPoints);
  boundaries->SetPolys(outputTriangles);

  boundaries->GetCellData()->AddArray(labelsBack);
  boundaries->GetCellData()->AddArray(labelsFront);
}

void regenerateBoundaries(vtkPoints *points, vtkFloatArray *labels,
			  vtkIntArray *connectivity, vtkShortArray *coords,
			  int currentTimeStep, std::vector<float3> &vertices,
			  std::vector<float3> &ivertices, std::vector<int> &indices,
			  std::vector<int> &splitTimes)
{
  float co[3][12];
  float coc[3][3];  
  float po[3][12];
  float poc[3][3];

  setupNodeOffsets(points, connectivity, co, coc, po, poc);

  // detect if reset of usedEdges is necessary by checking the number of vertices  
  static std::map<std::pair<int,int>, char> usedEdges;
  if (vertices.size() == 0) {
    usedEdges.clear();
  }

  int numPrevVertices = vertices.size();
  int vidx = 0;

  std::vector<float3> verticesTmp;
  std::vector<float3> iverticesTmp;

  const int numPoints = points->GetNumberOfPoints();
  for (int i = 0; i < numPoints; ++i) {

    int conn[3] = {connectivity->GetComponent(i, 0),
  		   connectivity->GetComponent(i, 1),
  		   connectivity->GetComponent(i, 2)};

    float l0 = labels->GetValue(i);
    float c0[3] = {coords->GetComponent(i, 0),
		   coords->GetComponent(i, 1),
		   coords->GetComponent(i, 2)};

    double p0[3];
    points->GetPoint(i, p0);
    
    for (int j = 0; j < 3; ++j) {
      if (conn[j] > -1) {

  	float l1 = labels->GetValue(conn[j]);
  	double p1[3];
  	points->GetPoint(conn[j], p1);
      
  	if (l0 != l1) {

	  if (usedEdges.find(std::pair<int,int>(conn[j],i)) == usedEdges.end()) {
	    usedEdges[std::pair<int,int>(conn[j],i)] = 1;
	      
	    splitTimes.push_back(currentTimeStep);
	    splitTimes.push_back(currentTimeStep);
	    splitTimes.push_back(currentTimeStep);
	    splitTimes.push_back(currentTimeStep);

	    float3 verts[5];
	    verts[4] = make_float3(p0[0]+poc[j][0], 
				   p0[1]+poc[j][1], 
				   p0[2]+poc[j][2]);	 
	    float3 iverts[5];
	    iverts[4] = make_float3(c0[0]+coc[j][0], 
				    c0[1]+coc[j][1], 
				    c0[2]+coc[j][2]);
		  
	  for (int k = 0; k < 4; ++k) {

	      float3 vertex = {p0[0]+po[j][k*3+0], 
			       p0[1]+po[j][k*3+1], 
			       p0[2]+po[j][k*3+2]};
	      verts[k] = vertex;

	      float3 ivertex = {c0[0]+co[j][k*3+0], 
				c0[1]+co[j][k*3+1], 
				c0[2]+co[j][k*3+2]};
	      iverts[k] = ivertex;
	    }
	    verticesTmp.push_back(verts[0]);
	    verticesTmp.push_back(verts[1]);
	    verticesTmp.push_back(verts[4]);
	    verticesTmp.push_back(verts[1]);
	    verticesTmp.push_back(verts[2]);
	    verticesTmp.push_back(verts[4]);
	    verticesTmp.push_back(verts[2]);
	    verticesTmp.push_back(verts[3]);
	    verticesTmp.push_back(verts[4]);
	    verticesTmp.push_back(verts[3]);
	    verticesTmp.push_back(verts[0]);
	    verticesTmp.push_back(verts[4]);
	    
	    iverticesTmp.push_back(iverts[0]);
	    iverticesTmp.push_back(iverts[1]);
	    iverticesTmp.push_back(iverts[4]);
	    iverticesTmp.push_back(iverts[1]);
	    iverticesTmp.push_back(iverts[2]);
	    iverticesTmp.push_back(iverts[4]);
	    iverticesTmp.push_back(iverts[2]);
	    iverticesTmp.push_back(iverts[3]);
	    iverticesTmp.push_back(iverts[4]);
	    iverticesTmp.push_back(iverts[3]);
	    iverticesTmp.push_back(iverts[0]);
	    iverticesTmp.push_back(iverts[4]);

	    vidx += 12;	  
	  }
	}
      }
    }
  }

  // TODO:
  // to generalize code, store the vertex map between function executions
  std::vector<int> vertexList;
  mergeVertices(iverticesTmp, vertices.size(), indices, vertexList);

  for (int i = 0; i < vertexList.size(); ++i) {
    int id = vertexList[i];    
    vertices.push_back(verticesTmp[id]);
    ivertices.push_back(iverticesTmp[id]);
  }
}

void regenerateBoundaries(vtkPoints *points, vtkFloatArray *labels,
			  vtkIntArray *connectivity, vtkShortArray *coords,
			  int currentTimeStep, vtkPolyData *boundaries)
{
  float co[3][12];
  float coc[3][3];  
  float po[3][12];
  float poc[3][3];

  setupNodeOffsets(points, connectivity, co, coc, po, poc);

  const int numVertices = boundaries->GetPoints()->GetNumberOfPoints();
  vtkFloatArray *splitTimes = 
    vtkFloatArray::SafeDownCast(boundaries->GetCellData()->GetArray("SplitTime"));

  // detect if reset of usedEdges is necessary by checking the number of vertices  
  static std::map<std::pair<int,int>, char> usedEdges;
  if (numVertices == 0) {
    usedEdges.clear();
  }

  int numPrevVertices = numVertices;
  int vidx = 0;

  std::vector<float3> verticesTmp;
  std::vector<float3> iverticesTmp;

  const int numPoints = points->GetNumberOfPoints();
  for (int i = 0; i < numPoints; ++i) {

    int conn[3] = {connectivity->GetComponent(i, 0),
  		   connectivity->GetComponent(i, 1),
  		   connectivity->GetComponent(i, 2)};

    float l0 = labels->GetValue(i);
    float c0[3] = {coords->GetComponent(i, 0),
		   coords->GetComponent(i, 1),
		   coords->GetComponent(i, 2)};

    double p0[3];
    points->GetPoint(i, p0);
    
    for (int j = 0; j < 3; ++j) {
      if (conn[j] > -1) {

  	float l1 = labels->GetValue(conn[j]);
  	double p1[3];
  	points->GetPoint(conn[j], p1);
      
  	if (l0 != l1) {

	  if (usedEdges.find(std::pair<int,int>(conn[j],i)) == usedEdges.end()) {
	    usedEdges[std::pair<int,int>(conn[j],i)] = 1;
	      
	    splitTimes->InsertNextTuple1(currentTimeStep);
	    splitTimes->InsertNextTuple1(currentTimeStep);
	    splitTimes->InsertNextTuple1(currentTimeStep);
	    splitTimes->InsertNextTuple1(currentTimeStep);

	    float3 verts[5];
	    verts[4] = make_float3(p0[0]+poc[j][0], 
				   p0[1]+poc[j][1], 
				   p0[2]+poc[j][2]);	 
	    float3 iverts[5];
	    iverts[4] = make_float3(c0[0]+coc[j][0], 
				    c0[1]+coc[j][1], 
				    c0[2]+coc[j][2]);
		  

	  for (int k = 0; k < 4; ++k) {

	      float3 vertex = {p0[0]+po[j][k*3+0], 
			       p0[1]+po[j][k*3+1], 
			       p0[2]+po[j][k*3+2]};
	      verts[k] = vertex;

	      float3 ivertex = {c0[0]+co[j][k*3+0], 
				c0[1]+co[j][k*3+1], 
				c0[2]+co[j][k*3+2]};
	      iverts[k] = ivertex;
	    }
	    verticesTmp.push_back(verts[0]);
	    verticesTmp.push_back(verts[1]);
	    verticesTmp.push_back(verts[4]);
	    verticesTmp.push_back(verts[1]);
	    verticesTmp.push_back(verts[2]);
	    verticesTmp.push_back(verts[4]);
	    verticesTmp.push_back(verts[2]);
	    verticesTmp.push_back(verts[3]);
	    verticesTmp.push_back(verts[4]);
	    verticesTmp.push_back(verts[3]);
	    verticesTmp.push_back(verts[0]);
	    verticesTmp.push_back(verts[4]);
	    
	    iverticesTmp.push_back(iverts[0]);
	    iverticesTmp.push_back(iverts[1]);
	    iverticesTmp.push_back(iverts[4]);
	    iverticesTmp.push_back(iverts[1]);
	    iverticesTmp.push_back(iverts[2]);
	    iverticesTmp.push_back(iverts[4]);
	    iverticesTmp.push_back(iverts[2]);
	    iverticesTmp.push_back(iverts[3]);
	    iverticesTmp.push_back(iverts[4]);
	    iverticesTmp.push_back(iverts[3]);
	    iverticesTmp.push_back(iverts[0]);
	    iverticesTmp.push_back(iverts[4]);
	    vidx += 12;	  
	  }
	}
      }
    }
  }

  // TODO:
  // to generalize code, store the vertex map between function executions
  std::vector<int> vertexList;
  vtkCellArray *cells = boundaries->GetPolys();
  

  mergeVertices(iverticesTmp, numVertices, cells, vertexList);
  
  vtkPoints *bpoints = boundaries->GetPoints();
  vtkFloatArray *ivertices = 
    vtkFloatArray::SafeDownCast(boundaries->GetPointData()->GetArray("IVertices"));

  for (int i = 0; i < vertexList.size(); ++i) {

    int id = vertexList[i];    

    bpoints->InsertNextPoint(verticesTmp[id].x, verticesTmp[id].y, verticesTmp[id].z);
    ivertices->InsertNextTuple3(iverticesTmp[id].x, iverticesTmp[id].y, iverticesTmp[id].z);
  }
}

void mergePatches(std::vector<float3> &vertices,
		  std::vector<float3> &ivertices,
		  std::vector<int> &indices)
{
  // find used vertices
  std::map<float3, int, compare_float3> vertexMap;
  vertexMap.clear();
  std::vector<int> useMask(ivertices.size());
  std::vector<float3> usedVertices(0);
  std::vector<int> newIndices(indices.size());


  for (int i = 0; i < ivertices.size(); ++i) {
    const float3 &key = ivertices[i];

    if (vertexMap.find(key) == vertexMap.end()) {	
      useMask[i] = usedVertices.size();
      vertexMap[key] = useMask[i];
      usedVertices.push_back(vertices[i]);
    }
    else {
      useMask[i] = vertexMap[key];
    }
  }

  for (int i = 0; i < indices.size(); ++i) {
    newIndices[i] = useMask[indices[i]];
  }
  vertices = usedVertices;
  indices = newIndices;
}

void mergePatches(vtkPolyData *boundaries)
{
  vtkFloatArray *ivertices =
    vtkFloatArray::SafeDownCast(boundaries->GetPointData()->GetArray("IVertices"));
  vtkPoints *vertices = boundaries->GetPoints();
  
  // find used vertices
  std::map<float3, int, compare_float3> vertexMap;
  vertexMap.clear();
  std::vector<int> useMask(ivertices->GetNumberOfTuples());
  vtkPoints *usedVertices = vtkPoints::New();
  
  for (int i = 0; i < ivertices->GetNumberOfTuples(); ++i) {
    const float3 &key = make_float3(ivertices->GetComponent(i,0),
				    ivertices->GetComponent(i,1),
				    ivertices->GetComponent(i,2));

    if (vertexMap.find(key) == vertexMap.end()) {	
      useMask[i] = usedVertices->GetNumberOfPoints();
      vertexMap[key] = useMask[i];
      double p[3];
      vertices->GetPoint(i, p);
      usedVertices->InsertNextPoint(p);
    }
    else {
      useMask[i] = vertexMap[key];
    }
  }

  vtkIdTypeArray *indices = boundaries->GetPolys()->GetData();
  int numTriangles = boundaries->GetPolys()->GetNumberOfCells();
  for (int t = 0; t < numTriangles; ++t) {
    indices->SetValue(t*4, 3);
    for (int v = 1; v <= 3; ++v) {
      indices->SetValue(t*4+v, useMask[indices->GetValue(t*4+v)]);//!!
    }
  }

  vertices->Delete();
  boundaries->SetPoints(usedVertices);
}
