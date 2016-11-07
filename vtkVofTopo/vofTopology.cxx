#include "vofTopology.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkShortArray.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLRectilinearGridWriter.h"

#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <cfloat>
#include <set>
#include <unordered_set>
#include <cmath>
#include <array>
#include <omp.h>
#include <algorithm>

#include "marchingCubes_cpu.h"

double g_emf0 = 0.000001;
double g_emf1 = 0.999999;

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

  void placeSeedsPLIC(vtkPoints *seeds,
		      const float cellCenter[3], 
		      const float cellSize[3],
		      const int refinement,
		      const int cellRes[3],
		      const float f,
		      const std::vector<float> &lstar,
		      const std::vector<float> &normalsInt,
		      const double bounds[6],
		      const int idx)
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
 	  }
	}
      }
    }    
  }

  void placeSeeds(vtkPoints *seeds,
		  const float cellCenter[3], 
		  const float cellSize[3],
		  const int refinement,
		  const int cellRes[3],
		  vtkRectilinearGrid *grid,
		  vtkDataArray *data,
		  const double bounds[6],
		  const int idx)
  {
    const float fcell = data->GetComponent(0,idx);
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

    for (int zr = 0; zr < subdiv; ++zr) {
      for (int yr = 0; yr < subdiv; ++yr) {
	for (int xr = 0; xr < subdiv; ++xr) {

	  float dx[3] = {originOffset[0] + xr*cellSizeTmp[0],
			 originOffset[1] + yr*cellSizeTmp[1],
			 originOffset[2] + zr*cellSizeTmp[2]};
	  float seed[3] = {cellCenter[0]+dx[0], // pos within cell
			   cellCenter[1]+dx[1],
			   cellCenter[2]+dx[2]};


	  double x[3] = {seed[0],seed[1],seed[2]};
	  double pcoords[3];
	  int ijk[3];
	  grid->ComputeStructuredCoordinates(x, ijk, pcoords);
	  
	  if (pointWithinBounds(seed, bounds)) {
	    if (fcell >= g_emf1) {
	      seeds->InsertNextPoint(seed);
	    }
	    else {
	      float f = interpolateScaCellBasedData(data, cellRes, ijk, pcoords);
	      if (f > g_emf0) {
		seeds->InsertNextPoint(seed);		  
	      }
	    }
 	  }
	}
      }
    }    
  }

  // data on nodes
  static float interpolateSca(const float *vofField,
			      const int* res, const int idxCell[3],
			      const double bcoords[3])
  {
    int lx = idxCell[0];
    int ly = idxCell[1];
    int lz = idxCell[2];
    float x = bcoords[0];
    float y = bcoords[1];
    float z = bcoords[2];
    
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
      vv[i] = vofField[id[i]];
    }

    float a = (1.0f-x)*vv[0] + x*vv[1];
    float b = (1.0f-x)*vv[2] + x*vv[3];
    float c = (1.0f-y)*a + y*b;
    a = (1.0f-x)*vv[4] + x*vv[5];
    b = (1.0f-x)*vv[6] + x*vv[7];
    float d = (1.0f-y)*a + y*b;

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

  int ComputeStructuredCoordinates(vtkRectilinearGrid *grid,
				 const double x[3], int ijk[3])
  {
    int i, j;
    double xPrev, xNext, tmp;
    vtkDataArray *scalars[3];
    int Dimensions[3];

    grid->GetDimensions(Dimensions);

    scalars[0] = grid->GetXCoordinates();
    scalars[1] = grid->GetYCoordinates();
    scalars[2] = grid->GetZCoordinates();
    //
    // Find locations in x-y-z direction
    //
    ijk[0] = ijk[1] = ijk[2] = 0;

    for ( j=0; j < 3; j++ )
      {
	xPrev = scalars[j]->GetComponent(0, 0);
	xNext = scalars[j]->GetComponent(scalars[j]->GetNumberOfTuples()-1, 0);
	if (xNext < xPrev)
	  {
	    tmp = xNext;
	    xNext = xPrev;
	    xPrev = tmp;
	  }
	if ( x[j] < xPrev || x[j] > xNext )
	  {
	    return 0;
	  }
	if (x[j] == xNext  && Dimensions[j] != 1)
	  {
	    return 0;
	  }

	for (i=1; i < scalars[j]->GetNumberOfTuples(); i++)
	  {
	    xNext = scalars[j]->GetComponent(i, 0);
	    if ( x[j] >= xPrev && x[j] < xNext )
	      {
		ijk[j] = i - 1;
		break;
	      }

	    else if ( x[j] == xNext )
	      {
		ijk[j] = i - 1;
		break;
	      }
	    xPrev = xNext;
	  }
      }

    return 1;
  }

}

float interpolateSca(vtkDataArray *vofField,
		     const int* res, const int idxCell[3],
		     const double bcoords[3])
{
  int lx = idxCell[0];
  int ly = idxCell[1];
  int lz = idxCell[2];
  float x = bcoords[0];
  float y = bcoords[1];
  float z = bcoords[2];
    
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
    vv[i] = vofField->GetComponent(id[i],0);
  }

  float a = (1.0f-x)*vv[0] + x*vv[1];
  float b = (1.0f-x)*vv[2] + x*vv[3];
  float c = (1.0f-y)*a + y*b;
  a = (1.0f-x)*vv[4] + x*vv[5];
  b = (1.0f-x)*vv[6] + x*vv[7];
  float d = (1.0f-y)*a + y*b;

  return (1.0f-z)*c + z*d;
}


void computeGradient(vtkDataArray *data,
		     const int cellRes[3], int ijk[3],
		     vtkFloatArray *coordCenters[3],
		     double pcoords[3], float grad[3])
{
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];
  int im = std::max(i-1,0);
  int ip = std::min(i+1,cellRes[0]-1);
  float dx = coordCenters[0]->GetComponent(ip,0) - coordCenters[0]->GetComponent(im,0);
  int jm = std::max(j-1,0);	  
  int jp = std::min(j+1,cellRes[1]-1);
  float dy = coordCenters[1]->GetComponent(jp,0) - coordCenters[1]->GetComponent(jm,0);
  int km = std::max(k-1,0);	  
  int kp = std::min(k+1,cellRes[2]-1);
  float dz = coordCenters[2]->GetComponent(kp,0) - coordCenters[2]->GetComponent(km,0);

  int id_left[3]   = {im, j, k};
  int id_right[3]  = {ip, j, k};
  int id_bottom[3] = {i, jm, k};
  int id_top[3]    = {i, jp, k};
  int id_back[3]   = {i, j, km};
  int id_front[3]  = {i, j, kp};

  float3 ddx = {dx,dy,dz};
  float3 dgrad = {
    (interpolateScaCellBasedData(data, cellRes, id_right, pcoords) - 
     interpolateScaCellBasedData(data, cellRes, id_left, pcoords))/dx,
    (interpolateScaCellBasedData(data, cellRes, id_top, pcoords) - 
     interpolateScaCellBasedData(data, cellRes, id_bottom, pcoords))/dy,
    (interpolateScaCellBasedData(data, cellRes, id_front, pcoords) - 
     interpolateScaCellBasedData(data, cellRes, id_back, pcoords))/dz};

  if (length(dgrad) > 0.0) {
    dgrad = normalize(dgrad);
  }
  dgrad *= ddx/2.0f;

  grad[0] = dgrad.x;
  grad[1] = dgrad.y;
  grad[2] = dgrad.z;
}

float interpolateScaCellBasedData(vtkDataArray *scalarField, const int* cellRes,
				  const int idxCell[3], const double bcoords[3])
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
  if (ux > cellRes[0]-1) ux = cellRes[0]-1;
  if (uy > cellRes[1]-1) uy = cellRes[1]-1;
  if (uz > cellRes[2]-1) uz = cellRes[2]-1;

  unsigned lzslab = lz*cellRes[0]*cellRes[1];
  unsigned uzslab = uz*cellRes[0]*cellRes[1];
  int lyr = ly*cellRes[0];
  int uyr = uy*cellRes[0];

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
    vv[i] = scalarField->GetComponent(id[i], 0);
  }

  float a = (1.0f-x)*vv[0] + x*vv[1];
  float b = (1.0f-x)*vv[2] + x*vv[3];
  float c = (1.0f-y)*a + y*b;
  a = (1.0f-x)*vv[4] + x*vv[5];
  b = (1.0f-x)*vv[6] + x*vv[7];
  float d = (1.0f-y)*a + y*b;

  return (1.0f-z)*c + z*d;
}

float3 interpolateVecCellBasedData(vtkDataArray *velocityField,
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

void computeL(const int cellRes[3],
	      const std::vector<float> &dx,
	      const std::vector<float> &dy,
	      const std::vector<float> &dz, 
	      vtkDataArray *f,
	      const std::vector<float> &normals,
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
	  if (i > 0 && f->GetComponent(fo-1,0) < g_emf0 || 
	      i < w-1 && f->GetComponent(fo+1,0) < g_emf0 || 
	      j > 0 && f->GetComponent(fo-w,0) < g_emf0 || 
	      j < h-1 && f->GetComponent(fo+w,0) < g_emf0 || 
	      k > 0 && f->GetComponent(fo-w*h,0) < g_emf0 || 
	      k < d-1 && f->GetComponent(fo+w*h,0) < g_emf0) {
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

void generateSeedPointsOnNodes(vtkRectilinearGrid *velocityGrid,
			       vtkPoints *points,
			       int globalExtent[6],
			       int numGhostLevels)
{
  vtkDataArray *coords[3] = {velocityGrid->GetXCoordinates(), 
			     velocityGrid->GetYCoordinates(), 
			     velocityGrid->GetZCoordinates()};

  int extent[6];
  velocityGrid->GetExtent(extent);

  int imin = extent[0] + (extent[0] > globalExtent[0] ? numGhostLevels : 0);
  int imax = extent[1] - (extent[1] < globalExtent[1] ? numGhostLevels : 0);
  int jmin = extent[2] + (extent[2] > globalExtent[2] ? numGhostLevels : 0);
  int jmax = extent[3] - (extent[3] < globalExtent[3] ? numGhostLevels : 0);
  int kmin = extent[4] + (extent[4] > globalExtent[4] ? numGhostLevels : 0);
  int kmax = extent[5] - (extent[5] < globalExtent[5] ? numGhostLevels : 0);

  for (int k = kmin; k <= kmax; ++k) {
    double z = coords[2]->GetComponent(k-extent[4],0);
    for (int j = jmin; j <= jmax; ++j) {
      double y = coords[1]->GetComponent(j-extent[2],0);
      for (int i = imin; i <= imax; ++i) {
	double x = coords[0]->GetComponent(i-extent[0],0);
	points->InsertNextPoint(x,y,z);
      }
    }
  }
}

void generateSeedPointsInCellCenters(vtkRectilinearGrid *velocityGrid,
				     vtkPoints *points,
				     int globalExtent[6],
				     int numGhostLevels)
{
  vtkDataArray *coords[3] = {velocityGrid->GetXCoordinates(), 
			     velocityGrid->GetYCoordinates(), 
			     velocityGrid->GetZCoordinates()};

  int extent[6];
  velocityGrid->GetExtent(extent);

  int imin = extent[0] + (extent[0] > globalExtent[0] ? numGhostLevels : 0);
  int imax = extent[1] - (extent[1] < globalExtent[1] ? numGhostLevels : 0);
  int jmin = extent[2] + (extent[2] > globalExtent[2] ? numGhostLevels : 0);
  int jmax = extent[3] - (extent[3] < globalExtent[3] ? numGhostLevels : 0);
  int kmin = extent[4] + (extent[4] > globalExtent[4] ? numGhostLevels : 0);
  int kmax = extent[5] - (extent[5] < globalExtent[5] ? numGhostLevels : 0);

  for (int k = kmin; k < kmax; ++k) {
    double z0 = coords[2]->GetComponent(k-extent[4],0);
    double z1 = coords[2]->GetComponent(k-extent[4]+1,0);
    double z = (z0+z1)/2.0;
    for (int j = jmin; j < jmax; ++j) {
      double y0 = coords[1]->GetComponent(j-extent[2],0);
      double y1 = coords[1]->GetComponent(j-extent[2]+1,0);
      double y = (y0+y1)/2.0;
      for (int i = imin; i < imax; ++i) {
	double x0 = coords[0]->GetComponent(i-extent[0],0);
	double x1 = coords[0]->GetComponent(i-extent[0]+1,0);
	double x = (x0+x1)/2.0;
	points->InsertNextPoint(x,y,z);
      }
    }
  }
}

void generateSeedPointsInCellCenters(vtkRectilinearGrid *velocityGrid,
				     int refinement,
				     std::vector<float4> &points,
				     int globalExtent[6],
				     int numGhostLevels)
{
  vtkDataArray *coords[3] = {velocityGrid->GetXCoordinates(), 
			     velocityGrid->GetYCoordinates(), 
			     velocityGrid->GetZCoordinates()};

  int extent[6];
  velocityGrid->GetExtent(extent);

  int imin = extent[0] + (extent[0] > globalExtent[0] ? numGhostLevels : 0);
  int imax = extent[1] - (extent[1] < globalExtent[1] ? numGhostLevels : 0);
  int jmin = extent[2] + (extent[2] > globalExtent[2] ? numGhostLevels : 0);
  int jmax = extent[3] - (extent[3] < globalExtent[3] ? numGhostLevels : 0);
  int kmin = extent[4] + (extent[4] > globalExtent[4] ? numGhostLevels : 0);
  int kmax = extent[5] - (extent[5] < globalExtent[5] ? numGhostLevels : 0);

  for (int k = kmin; k <= kmax; ++k) {
    double z0 = coords[2]->GetComponent(k-extent[4],0);
    double z1 = coords[2]->GetComponent(k-extent[4]+1,0);
    double z = (z0+z1)/2.0;
    for (int j = jmin; j <= jmax; ++j) {
      double y0 = coords[1]->GetComponent(j-extent[2],0);
      double y1 = coords[1]->GetComponent(j-extent[2]+1,0);
      double y = (y0+y1)/2.0;
      for (int i = imin; i <= imax; ++i) {
	double x0 = coords[0]->GetComponent(i-extent[0],0);
	double x1 = coords[0]->GetComponent(i-extent[0]+1,0);
	double x = (x0+x1)/2.0;
	points.push_back(make_float4(x,y,z,1.0f));
      }
    }
  }
}


void generateSeedPointsInCellCenters(vtkRectilinearGrid *velocityGrid,
				     int refinement,
				     std::vector<float4> &points,
				     int particleExtent[6],
				     int globalExtent[6],
				     int numGhostLevels)
{
  vtkDataArray *coords[3] = {velocityGrid->GetXCoordinates(), 
			     velocityGrid->GetYCoordinates(), 
			     velocityGrid->GetZCoordinates()};

  int imin = particleExtent[0];
  int imax = particleExtent[1];
  int jmin = particleExtent[2];
  int jmax = particleExtent[3];
  int kmin = particleExtent[4];
  int kmax = particleExtent[5];

  for (int k = kmin; k <= kmax; ++k) {
    double z0 = coords[2]->GetComponent(k,0);
    double z1 = coords[2]->GetComponent(k+1,0);
    double z = (z0+z1)/2.0;
    for (int j = jmin; j <= jmax; ++j) {
      double y0 = coords[1]->GetComponent(j,0);
      double y1 = coords[1]->GetComponent(j+1,0);
      double y = (y0+y1)/2.0;
      for (int i = imin; i <= imax; ++i) {
	double x0 = coords[0]->GetComponent(i,0);
	double x1 = coords[0]->GetComponent(i+1,0);
	double x = (x0+x1)/2.0;
	points.push_back(make_float4(x,y,z,1.0f));
      }
    }
  }
}

float neighborF(vtkDataArray *data, int i, int j, int k, const int cellRes[3])
{
  int im = std::max(i-1,0);
  int ip = std::min(i+1,cellRes[0]-1);
  int jm = std::max(j-1,0);	  
  int jp = std::min(j+1,cellRes[1]-1);
  int km = std::max(k-1,0);	  
  int kp = std::min(k+1,cellRes[2]-1);

  float ff = 0.0f;
  int idx;

  idx = im + j*cellRes[0] + k*cellRes[0]*cellRes[1];
  ff += data->GetComponent(idx,0);
  idx = ip + j*cellRes[0] + k*cellRes[0]*cellRes[1];
  ff += data->GetComponent(idx,0);
  idx = i + jm*cellRes[0] + k*cellRes[0]*cellRes[1];
  ff += data->GetComponent(idx,0);
  idx = i + jp*cellRes[0] + k*cellRes[0]*cellRes[1];
  ff += data->GetComponent(idx,0);
  idx = i + j*cellRes[0] + km*cellRes[0]*cellRes[1];
  ff += data->GetComponent(idx,0);
  idx = i + j*cellRes[0] + kp*cellRes[0]*cellRes[1];
  ff += data->GetComponent(idx,0);

  return ff;
}

void generateSeedPoints(vtkRectilinearGrid *vofGrid,
			int refinement,
			vtkPoints *points,
			int globalExtent[6],
			int numGhostLevels,
			int plicSeeding)
{
  vtkDataArray *vofArray = vofGrid->GetCellData()->GetArray(0);
  
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
  if (plicSeeding) {
    normals.resize(nodeRes[0]*nodeRes[1]*nodeRes[2]*3);
    computeNormals(nodeRes, dx[0], dx[1], dx[2], vofArray, normals);
  }
  std::vector<float> lstar;
  std::vector<float> normalsInt;
  if (plicSeeding) {
    lstar.resize(cellRes[0]*cellRes[1]*cellRes[2]);
    normalsInt.resize(cellRes[0]*cellRes[1]*cellRes[2]*3);
    computeL(cellRes, dx[0], dx[1], dx[2], vofArray, normals, lstar, normalsInt);
  }
  //--------------
  vtkDataArray *data = vofGrid->GetCellData()->GetArray(0);
  
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

	if (plicSeeding) {
	  if (f > g_emf0) {
	    placeSeedsPLIC(points, cellCenter, cellSize, refinement, cellRes,
			   f, lstar, normalsInt, bounds, idx);
	  }
	}
	// else {	  
	//   if (f > g_emf0) {
	//     placeSeeds(points, cellCenter, cellSize, refinement, cellRes, vofGrid, data, bounds, idx);
	//   }
	//   else {
	//     float fn = neighborF(data,i,j,k,cellRes);
	//     if (fn > g_emf0) {
	//       placeSeeds(points, cellCenter, cellSize, refinement, cellRes, vofGrid, data, bounds, idx);
	//     }
	//   }
	// }

	++icur;
      }
      ++jcur;
    }
    ++kcur;
  }
}


float4 vofCorrector(const float4 pos1, vtkDataArray *vofField,
		    vtkDataArray *coords[3], const int res[3],
		    const int ijk[3],
		    vtkRectilinearGrid *vofGrid, float &fCorrected)
{
  int stencilRange = 4;  
  int x = ijk[0];
  int y = ijk[1];
  int z = ijk[2];

  int width = 2*stencilRange+1;
  int neighCnt = width*width*width;

  double3 xyzi = {
    (coords[0]->GetComponent(x+1,0) + coords[0]->GetComponent(x,0))/2.0,
    (coords[1]->GetComponent(y+1,0) + coords[1]->GetComponent(y,0))/2.0,
    (coords[2]->GetComponent(z+1,0) + coords[2]->GetComponent(z,0))/2.0};

  double mindist2 = std::numeric_limits<float>::max();
  int bestIdx[3] = {-1};
    
  for (int k = 0; k < neighCnt; k++) {
    
    int offsetX = (k%width) - stencilRange;
    int offsetY = ((k/width)%width) - stencilRange;
    int offsetZ = (k/(width*width)) - stencilRange;
    int neighborX = x + offsetX;
    int neighborY = y + offsetY;
    int neighborZ = z + offsetZ;
    
    if (neighborX >= 0 && neighborX < res[0] &&
	neighborY >= 0 && neighborY < res[1] &&
	neighborZ >= 0 && neighborZ < res[2])	{

      double3 xyzj = {
	(coords[0]->GetComponent(neighborX+1,0) + coords[0]->GetComponent(neighborX,0))/2.0,
	(coords[1]->GetComponent(neighborY+1,0) + coords[1]->GetComponent(neighborY,0))/2.0,
	(coords[2]->GetComponent(neighborZ+1,0) + coords[2]->GetComponent(neighborZ,0))/2.0};
      
      double X = xyzj.x - xyzi.x;
      double Y = xyzj.y - xyzi.y;
      double Z = xyzj.z - xyzi.z;

      int neighborIndex = neighborX + neighborY*res[0] + neighborZ*res[0]*res[1];
      double f = vofField->GetComponent(neighborIndex, 0);
      double dist2 = X*X + Y*Y + Z*Z;
      
      if (f > g_emf0 && dist2 <= mindist2) {
	bestIdx[0] = neighborX;
	bestIdx[1] = neighborY;
	bestIdx[2] = neighborZ;
	fCorrected = f;
	mindist2 = dist2;
      }      
    }
  }

  float4 pos2 = pos1;
  if (bestIdx[0] > -1) {
    pos2 = make_float4((coords[0]->GetComponent(bestIdx[0]+1,0) + coords[0]->GetComponent(bestIdx[0],0))/2.0,
		       (coords[1]->GetComponent(bestIdx[1]+1,0) + coords[1]->GetComponent(bestIdx[1],0))/2.0,
		       (coords[2]->GetComponent(bestIdx[2]+1,0) + coords[2]->GetComponent(bestIdx[2],0))/2.0,
		       1.0f);
  }
  return pos2;
}

float getCellVof(const double pos[3], vtkRectilinearGrid *grid,
		 vtkDataArray *vofArray1, const int cellRes[3],
		 int ijk[3])
{
  ComputeStructuredCoordinates(grid, pos, ijk);
  int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
  return vofArray1->GetComponent(idx, 0);
}




const static float3 ones = {1.0f,1.0f,1.0f};

float3 intersectBox(float3 pos, float3 dir, float3 boxmin, float3 boxmax)
{
  // compute intersection of ray with all six bbox planes
  float3 invR = ones / dir;
  float3 tbot = invR * (boxmin - pos);
  float3 ttop = invR * (boxmax - pos);

  // re-order intersections to find smallest and largest on each axis
  float3 tmin = minf3(ttop, tbot);
  float3 tmax = maxf3(ttop, tbot);
  float tnear = std::max(std::max(tmin.x, tmin.y), std::max(tmin.x, tmin.z));
  float w = std::min(std::min(boxmax.x-boxmin.x, boxmax.y-boxmin.y), boxmax.z-boxmin.z);
  tnear += w/16.0f;

  float3 pp = pos + dir*tnear;

  return pos + dir*tnear;
}

float4 vofCorrector2(const float4 pos1, vtkDataArray *vofField,
		     vtkDataArray *coords[3], const int res[3],
		     int ijk[3],
		     vtkRectilinearGrid *vofGrid, float &fCorrected)
{
  int stencilRange = 4;
  int x = ijk[0];
  int y = ijk[1];
  int z = ijk[2];

  int width = 2*stencilRange+1;
  int neighCnt = width*width*width;

  float3 xyzi = {
    (coords[0]->GetComponent(x+1,0) + coords[0]->GetComponent(x,0))/2.0,
    (coords[1]->GetComponent(y+1,0) + coords[1]->GetComponent(y,0))/2.0,
    (coords[2]->GetComponent(z+1,0) + coords[2]->GetComponent(z,0))/2.0};

  double mindist2 = std::numeric_limits<float>::max();
  int bestIdx[3] = {-1};

  for (int k = 0; k < neighCnt; k++) {
    
    int offsetX = (k%width) - stencilRange;
    int offsetY = ((k/width)%width) - stencilRange;
    int offsetZ = (k/(width*width)) - stencilRange;
    int neighborX = x + offsetX;
    int neighborY = y + offsetY;
    int neighborZ = z + offsetZ;
    
    if (neighborX >= 0 && neighborX < res[0] &&
	neighborY >= 0 && neighborY < res[1] &&
	neighborZ >= 0 && neighborZ < res[2])	{

      float3 xyzj = {
	(coords[0]->GetComponent(neighborX+1,0) + coords[0]->GetComponent(neighborX,0))/2.0,
	(coords[1]->GetComponent(neighborY+1,0) + coords[1]->GetComponent(neighborY,0))/2.0,
	(coords[2]->GetComponent(neighborZ+1,0) + coords[2]->GetComponent(neighborZ,0))/2.0};
      
      double X = xyzj.x - xyzi.x;
      double Y = xyzj.y - xyzi.y;
      double Z = xyzj.z - xyzi.z;

      int neighborIndex = neighborX + neighborY*res[0] + neighborZ*res[0]*res[1];
      double f = vofField->GetComponent(neighborIndex, 0);
      double dist2 = X*X + Y*Y + Z*Z;
      
      if (f > g_emf0 && dist2 <= mindist2) {
	bestIdx[0] = neighborX;
	bestIdx[1] = neighborY;
	bestIdx[2] = neighborZ;
	fCorrected = f;
	mindist2 = dist2;
      }      
    }
  }

  float4 pos2 = pos1;  
  if  (bestIdx[0] > -1) {
    float3 boxmin = make_float3(coords[0]->GetComponent(bestIdx[0],0),
    				coords[1]->GetComponent(bestIdx[1],0),
    				coords[2]->GetComponent(bestIdx[2],0));
    float3 boxmax = make_float3(coords[0]->GetComponent(bestIdx[0]+1,0),
    				coords[1]->GetComponent(bestIdx[1]+1,0),
    				coords[2]->GetComponent(bestIdx[2]+1,0));
    float3 cellCenter = (boxmin + boxmax)/2.0f;
    float3 pos3f = make_float3(pos1);
    float3 pos4 = intersectBox(pos3f, normalize(cellCenter-pos3f), boxmin, boxmax);

    pos2 = make_float4(pos4);
    pos2.w = fCorrected;
    ijk[0] = bestIdx[0];
    ijk[1] = bestIdx[1];
    ijk[2] = bestIdx[2];

  }
  return pos2;
}

// float4 plicCorrector(const float4 pos1, const float3 normal,
// 		     const float lstar, float coords[6])
// {
//   float4 pos2 = pos1;
//   float3 nodes[8] = {make_float3(coords[0],coords[2],coords[4]),
// 		     make_float3(coords[1],coords[2],coords[4]),
// 		     make_float3(coords[0],coords[3],coords[4]),
// 		     make_float3(coords[1],coords[3],coords[4]),
// 		     make_float3(coords[0],coords[2],coords[5]),
// 		     make_float3(coords[1],coords[2],coords[5]),
// 		     make_float3(coords[0],coords[3],coords[5]),
// 		     make_float3(coords[1],coords[3],coords[5])};   
//   float3 attachPoint =
//     make_float3(normal.x > 0 ? coords[0] : coords[1],
// 		normal.y > 0 ? coords[2] : coords[3],
// 		normal.z > 0 ? coords[4] : coords[5]);
//   float projNodes[8];
//   for (int i = 0; i < 8; ++i) {
//     float3 posVec = nodes[i] - attachPoint;
//     projNodes[i] = std::abs(dot(posVec,normal));
//   }
//   char edges[12][2] = {{0,1},{0,2},{1,3},{2,3},
// 		       {4,5},{4,6},{5,7},{6,7},
// 		       {0,4},{1,5},{2,6},{3,7}};
  
//   std::vector<float3> edgeisects;
//   for (int i = 0; i < 12; ++i) {
//     int e0 = edges[i][0];
//     int e1 = edges[i][1];
//     float f0 = projNodes[e0];
//     float f1 = projNodes[e1];
//     if (f0 == f1) continue;
//     float a = (lstar-f0)/(f1-f0);
//     if (a > 0 && a < 1) {
//       edgeisects.push_back(nodes[e0]*(1.0f-a) + nodes[e1]*a);
//     }
//   }

//   float3 posVec = make_float3(pos1) - attachPoint;
//   float d = std::abs(dot(posVec,normal));

//   if (d > lstar) {
//     float d2 = d - lstar;
//     float4 disp = make_float4(d2*normal);
//     pos2 -= disp;
//     if (pos2.x < coords[0] || pos2.x > coords[1] ||
// 	pos2.y < coords[2] || pos2.y > coords[3] ||
// 	pos2.z < coords[4] || pos2.z > coords[5]) {
//       float mindist = 10000.0f;
//       int isel = -1;
//       for (int i = 0; i < edgeisects.size(); ++i) {
// 	float dist = length(make_float3(pos1)-edgeisects[i]);
// 	if (dist < mindist) {
// 	  mindist = dist;
// 	  isel = i;
// 	}
//       }
//       pos2 = make_float4(edgeisects[isel]);
//       pos2.w = 1.0f;
//     }
//   }
//   return pos2;
// }

float4 plicCorrector(const float4 pos1, const float3 normal,
		     const float lstar, float coords[6])
{
  float4 pos2 = pos1;
  float3 attachPoint =
    make_float3(normal.x > 0 ? coords[0] : coords[1],
		normal.y > 0 ? coords[2] : coords[3],
		normal.z > 0 ? coords[4] : coords[5]);

  float3 posVec = make_float3(pos1) - attachPoint;
  float d = std::abs(dot(posVec,normal));

  if (d > lstar) {
    float3 posVecN = make_float3(0.0f);
    if (length(posVec) > 0.0f) {
      posVecN = normalize(posVec);
    }
    float ca = std::abs(dot(posVecN, normal));
    float e = (d-lstar)/ca;    
    float4 disp = make_float4(e*posVecN,0.0f);
    pos2 -= disp;
  }
  return pos2;
}


void correctParticles(std::vector<float4> &particles,
		      std::vector<float> &uncertainty,
		      vtkRectilinearGrid *vofGrid[2], vtkDataArray *vofArray1,
		      vtkDataArray *coords[3], const int cellRes[3],
		      int plicCorrection, int vofCorrection)
{
  if (plicCorrection == 0 && vofCorrection == 0) {
    return;
  }
  
  std::vector<std::vector<float>> dx(3);
  dx[0].resize(cellRes[0]);
  dx[1].resize(cellRes[1]);
  dx[2].resize(cellRes[2]);
  for (int c = 0; c < 3; ++c) {
    for (int i = 0; i < cellRes[c]; ++i) {
      dx[c][i] = coords[c]->GetComponent(i+1,0) - coords[c]->GetComponent(i,0);
    }
  }

  int nodeRes[3] = {cellRes[0]+1,cellRes[1]+1,cellRes[2]+1};
  
  std::vector<float> normalsNodes;
  std::vector<float> lstar;
  std::vector<float> normals;

  if (plicCorrection) {

    normalsNodes.resize(nodeRes[0]*nodeRes[1]*nodeRes[2]*3);
    lstar.resize(cellRes[0]*cellRes[1]*cellRes[2]);
    normals.resize(cellRes[0]*cellRes[1]*cellRes[2]*3);
    
    computeNormals(nodeRes, dx[0], dx[1], dx[2], vofArray1, normalsNodes);
    computeL(cellRes, dx[0], dx[1], dx[2], vofArray1, normalsNodes, lstar, normals);
  }

#pragma omp parallel for
  for (int i = 0; i < particles.size(); ++i) {

    if (particles[i].w > -1.0f) {

      int ijk[3];
      double p[3] = {particles[i].x, particles[i].y, particles[i].z};
      float f = getCellVof(p, vofGrid[1], vofArray1, cellRes, ijk);
      
      if (vofCorrection && f <= g_emf0) {
	float4 prev_pos1 = particles[i];
	particles[i] = vofCorrector(particles[i], vofArray1, coords, cellRes, ijk, vofGrid[1], f);
	uncertainty[i] += length(make_float3(particles[i] - prev_pos1));
      }
      if (plicCorrection && (f > g_emf0 && f < g_emf1)) {
      	double x[3] = {particles[i].x, particles[i].y, particles[i].z};
      	double pcoords[3];
      	vofGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
      	int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      	float cubeCoords[6] = {coords[0]->GetComponent(ijk[0],0),
      			       coords[0]->GetComponent(ijk[0]+1,0),
      			       coords[1]->GetComponent(ijk[1],0),
      			       coords[1]->GetComponent(ijk[1]+1,0),
      			       coords[2]->GetComponent(ijk[2],0),
      			       coords[2]->GetComponent(ijk[2]+1,0)};
      	float3 norm = make_float3(normals[idx*3+0],
      				  normals[idx*3+1],
      				  normals[idx*3+2]);
      	float4 prev_pos1 = particles[i];

      	particles[i] = plicCorrector(particles[i], norm, lstar[idx], cubeCoords);
      	uncertainty[i] += length(make_float3(particles[i] - prev_pos1));
      }
    }
  }
}

//==============================================================================
class IntKeys3 {
public:
  IntKeys3(short k1, short k2, short k3) : key1(k1), key2(k2), key3(k3) { }

  bool operator<(const IntKeys3 &right) const {
    return (key1 < right.key1 ||
            (key1 == right.key1 &&
             (key2 < right.key2 ||
              (key2 == right.key2 &&
               (key3 < right.key3))))
            );
  }

  short key1, key2, key3;
};

void mapCellsToParticles(vtkRectilinearGrid *grid,
			 std::vector<float4> &particles,
			 std::map<IntKeys3,std::vector<int>> &particlesInCells)
{
  int pIdx = 0;
  for (const auto &particle : particles) {
    
    double p[3] = {particle.x, particle.y, particle.z};
    int ijk[3];
    ComputeStructuredCoordinates(grid, p, ijk);
    
    IntKeys3 key(ijk[0],ijk[1],ijk[2]);
    
    if (particlesInCells.find(key) == particlesInCells.end()) {
      particlesInCells[key] = {pIdx};
    }
    else {
      particlesInCells[key].push_back(pIdx);
    }
    
    ++pIdx;
  }

}

void correctParticles2(std::vector<float4> &particles,
		       std::vector<float4> &oldParticles,
		       std::vector<float> &uncertainty,
		       vtkRectilinearGrid *vofGrid[2], vtkDataArray *vofArray1,
		       vtkDataArray *coords[3], const int cellRes[3],
		       int plicCorrection, int vofCorrection, int smartCorrection)
{
  if (plicCorrection == 0 && vofCorrection == 0) {
    return;
  }

  vtkRectilinearGrid *grid_t1 = vofGrid[1];
  
  std::vector<std::vector<float>> dx(3);
  if (plicCorrection) {
    dx[0].resize(cellRes[0]);
    dx[1].resize(cellRes[1]);
    dx[2].resize(cellRes[2]);
    for (int c = 0; c < 3; ++c) {
      for (int i = 0; i < cellRes[c]; ++i) {
	dx[c][i] = coords[c]->GetComponent(i+1,0) - coords[c]->GetComponent(i,0);
      }
    }
  }
  
  int nodeRes[3] = {cellRes[0]+1,cellRes[1]+1,cellRes[2]+1};
  
  std::vector<float> normalsNodes;
  std::vector<float> lstar;
  std::vector<float> normals;

  if (plicCorrection) {

    normalsNodes.resize(nodeRes[0]*nodeRes[1]*nodeRes[2]*3);
    lstar.resize(cellRes[0]*cellRes[1]*cellRes[2]);
    normals.resize(cellRes[0]*cellRes[1]*cellRes[2]*3);
    
    computeNormals(nodeRes, dx[0], dx[1], dx[2], vofArray1, normalsNodes);
    computeL(cellRes, dx[0], dx[1], dx[2], vofArray1, normalsNodes, lstar, normals);
  }

  std::vector<std::pair<int,float4>> updatedParticles;
  updatedParticles.clear();

  std::map<IntKeys3,std::vector<int>> particlesInCells;
  particlesInCells.clear();
  
  if (smartCorrection) {
    mapCellsToParticles(grid_t1, oldParticles, particlesInCells);
  }  

  const int cellRange = 1;

  // go over new positions -----------------------------------------
  int pIdx = 0;
  for (const auto &particle : particles) {  

    float4 updatedParticle = particle;
    
    int ijk[3];
    double x[3] = {updatedParticle.x, updatedParticle.y, updatedParticle.z};      
    const float f = getCellVof(x, grid_t1, vofArray1, cellRes, ijk);    
    const bool pushParticle = f < g_emf1;

    if (smartCorrection && f <= g_emf0) {

      double p[3] = {oldParticles[pIdx].x, oldParticles[pIdx].y, oldParticles[pIdx].z};      
      ComputeStructuredCoordinates(grid_t1, p, ijk);
	
      float minDist = std::numeric_limits<float>::max();
      int bestId = -1;
	
      for (int z = ijk[2]-cellRange; z <= ijk[2]+cellRange; ++z) {
    	for (int y = ijk[2]-cellRange; y <= ijk[2]+cellRange; ++y) {
    	  for (int x = ijk[0]-cellRange; x <= ijk[0]+cellRange; ++x) {

    	    IntKeys3 key(x,y,z);
    	    const std::vector<int> &nIds = particlesInCells[key];
	      
    	    for (const auto &nId : nIds) {
    	      if (nId != pIdx) {
    		int ijk2[3];
    		double p2[3] = {particles[nId].x, particles[nId].y, particles[nId].z};
    		float f2 = getCellVof(p2, grid_t1, vofArray1, cellRes, ijk2);
    		if (f2 > g_emf0) {
		  
    		  float4 d = oldParticles[pIdx] - oldParticles[nId];
    		  float dist = length(make_float3(d.x,d.y,d.z));
    		  if (minDist > dist) {
    		    minDist = dist;
    		    bestId = nId;
    		  }
    		}
    	      }
    	    }
    	  }
    	}
      }		  
      if (bestId > -1) {
    	updatedParticle += particles[bestId] - oldParticles[bestId];
      }
    }
    if (pushParticle) {
      updatedParticles.push_back(std::pair<int,float4>(pIdx,updatedParticle));
    }
    ++pIdx;
  }

#pragma omp parallel for
  for (int i = 0; i < updatedParticles.size(); ++i) {

    float4 &updatedParticle = updatedParticles[i].second;

    double p[3] = {updatedParticle.x, updatedParticle.y, updatedParticle.z};
    int ijk[3];
    float f = getCellVof(p, grid_t1, vofArray1, cellRes, ijk);

    if (vofCorrection && f <= g_emf0) {
      updatedParticle = vofCorrector2(updatedParticle, vofArray1, coords, cellRes, ijk, grid_t1, f);
    }
    if (plicCorrection && (f > g_emf0 && f < g_emf1)) {

      int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      float cubeCoords[6] = {coords[0]->GetComponent(ijk[0],0),
    			     coords[0]->GetComponent(ijk[0]+1,0),
    			     coords[1]->GetComponent(ijk[1],0),
    			     coords[1]->GetComponent(ijk[1]+1,0),
    			     coords[2]->GetComponent(ijk[2],0),
    			     coords[2]->GetComponent(ijk[2]+1,0)};
      float3 norm = make_float3(normals[idx*3+0],
    				normals[idx*3+1],
    				normals[idx*3+2]);
      
      updatedParticle = plicCorrector(updatedParticle, norm, lstar[idx], cubeCoords);
    }
  }

  for (int i = 0; i < updatedParticles.size(); ++i) {
    int pIdx = updatedParticles[i].first;
    uncertainty[pIdx] += length(make_float3(particles[pIdx] - updatedParticles[i].second));
    particles[pIdx] = updatedParticles[i].second;
  }
}

const float inv6 = 0.16666666666666666f;

float4 rungeKutta4Int(const float4 &pos0, vtkRectilinearGrid *velocityGrid[2],
		      vtkDataArray *velocityArray0, vtkDataArray *velocityArray1,
		      vtkDataArray *velocityArray2, const int cellRes[3], const float deltaT)
{
  float4 pos = pos0;
  float dt = deltaT;
  double x[3] = {pos.x, pos.y, pos.z};
  int ijk[3];
  double pcoords[3];
  velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);

  float4 k1 = make_float4(interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords),0.0f);
  float4 k1p = pos + k1*dt*0.5f;
  
  x[0] = k1p.x; x[1] = k1p.y; x[2] = k1p.z;
  velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
  float4 k2 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);
  float4 k2p = pos + k2*dt*0.5f;
  
  x[0] = k2p.x; x[1] = k2p.y; x[2] = k2p.z;
  velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
  float4 k3 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);
  float4 k3p = pos + k3*dt;

  x[0] = k3p.x; x[1] = k3p.y; x[2] = k3p.z;
  velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);    
  float4 k4 = make_float4(interpolateVecCellBasedData(velocityArray2, cellRes, ijk, pcoords),0.0f);

  pos += dt*inv6*(k1 + 2.0f*k2 + 2.0f*k3 + k4);

  return pos;
}

float4 rungeKutta4(const float4 &pos0, vtkRectilinearGrid *velocityGrid[2],
		   vtkDataArray *velocityArray0, vtkDataArray *velocityArray1,
		   const int cellRes[3], const float t, const float incr,
		   const float t0, const float t1, int numSubSteps)
{
  const int numSteps = numSubSteps;
  float4 pos = pos0;
  float ct = t;
  
  float DeltaT = t1 - t0; // > 0
  float et = incr > 0 ? t1 : t0;
  float deltaT = et - ct; // <> 0
  float ctr = ct - t0; // > 0
  float a;
  
  for (int i = 0; i < numSteps; ++i) {

    float dt = i < numSteps-1 ? deltaT/numSteps : et - ct;
    double x[3] = {pos.x, pos.y, pos.z};
    int ijk[3];
    double pcoords[3];
    velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);

    a = ctr/DeltaT;
    float4 k10 = make_float4(interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords),0.0f);
    float4 k11 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);    
    float4 k1 = lerp(k10,k11,a);
    float4 k1p = pos + k1*dt*0.5f;
    x[0] = k1p.x; x[1] = k1p.y; x[2] = k1p.z;
    velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);

    a = (ctr+dt*0.5f)/DeltaT;
    float4 k20 = make_float4(interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords),0.0f);
    float4 k21 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);
    float4 k2 = lerp(k20,k21,a);
    float4 k2p = pos + k2*dt*0.5f;
    x[0] = k2p.x; x[1] = k2p.y; x[2] = k2p.z;    
    velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);

    a = (ctr+dt*0.5f)/DeltaT;
    float4 k30 = make_float4(interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords),0.0f);
    float4 k31 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);
    float4 k3 = lerp(k30,k31,a);
    float4 k3p = pos + k3*dt;
    x[0] = k3p.x; x[1] = k3p.y; x[2] = k3p.z;
    velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);

    a = (ctr+dt)/DeltaT;
    float4 k40 = make_float4(interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords),0.0f);
    float4 k41 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);
    float4 k4 = lerp(k40,k41,a);

    pos += dt*inv6*(k1 + 2.0f*k2 + 2.0f*k3 + k4);
    ct += dt;
    ctr += dt; 
  }
  return pos;
}

float4 iterativeHeun(const float4 &pos0, vtkRectilinearGrid *velocityGrid[2],
		     vtkDataArray *velocityArray0, vtkDataArray *velocityArray1,
		     vtkDataArray *vofArray1, const int cellRes[3], const float deltaT)
{
  double x[3] = {pos0.x, pos0.y, pos0.z};
  int ijk[3];
  double pcoords[3];
  velocityGrid[0]->ComputeStructuredCoordinates(x, ijk, pcoords);
  float4 velocity0 = make_float4(interpolateVecCellBasedData(velocityArray0, cellRes, ijk, pcoords),0.0f);
  float4 pos1 = pos0 + deltaT*velocity0;     // initial guess - forward Euler
  float4 velocity1;

  const int maxNumIter = 20;
  for (int i = 0; i < maxNumIter; ++i) {

    x[0] = pos1.x;
    x[1] = pos1.y;
    x[2] = pos1.z;
      
    velocityGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);
    velocity1 = make_float4(interpolateVecCellBasedData(velocityArray1, cellRes, ijk, pcoords),0.0f);

    float4 velocity = (velocity0 + velocity1)/2.0f;
    pos1 = pos0 + deltaT*velocity;	
  }      

  x[0] = pos1.x;
  x[1] = pos1.y;
  x[2] = pos1.z;
  int particleInsideGrid = velocityGrid[1]->ComputeStructuredCoordinates(x, ijk, pcoords);        
  if (particleInsideGrid) {
    int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
    pos1.w = vofArray1->GetComponent(idx, 0);
  }
  return pos1;
}

// iterative, solved with fixed point method - Newton's method can be viewed as such
// https://en.wikipedia.org/wiki/Fixed-point_iteration
// https://en.wikipedia.org/wiki/Trapezoidal_rule_%28differential_equations%29
void advectParticles(vtkRectilinearGrid *vofGrid[2],
		     vtkRectilinearGrid *velocityGrid[2],
		     std::vector<float4> &particles,
		     std::vector<float> &uncertainty,
		     const float deltaT, int integrationMethod,
		     int plicCorrection, int vofCorrection,
		     int smartCorrection, int numSubSteps)
{
  vtkDataArray *velocityArray0 = velocityGrid[0]->GetCellData()->GetArray(0);
  vtkDataArray *velocityArray1 = velocityGrid[1]->GetCellData()->GetArray(0);
  vtkDataArray *vofArray1 = vofGrid[1]->GetCellData()->GetArray(0);

  int nodeRes[3];
  vofGrid[1]->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  vtkDataArray *coords[3] = {vofGrid[1]->GetXCoordinates(),
			     vofGrid[1]->GetYCoordinates(),
			     vofGrid[1]->GetZCoordinates()};

  std::vector<float4> oldParticles = particles;
  
#pragma omp parallel for  
  for (int i = 0; i < particles.size(); ++i) {

    // if (particles[i].w <= g_emf0) {
    //   continue;
    // }
    if (integrationMethod == 0) { // Heun
      particles[i] = iterativeHeun(particles[i], velocityGrid, velocityArray0, velocityArray1,
				   vofArray1, cellRes, deltaT);
    }
    else {
      particles[i] = rungeKutta4(particles[i], velocityGrid, velocityArray0, velocityArray1, cellRes,
				 0.0f, 1.0f, 0.0f, deltaT, numSubSteps);
    }
  }

  //correctParticles(particles, uncertainty, vofGrid, vofArray1, coords, cellRes, plicCorrection, vofCorrection);
  correctParticles2(particles, oldParticles, uncertainty, vofGrid, vofArray1, coords, cellRes, plicCorrection, vofCorrection, smartCorrection);
}

// void advectParticles(vtkRectilinearGrid *velocityGrid[2],
// 		     std::vector<float4> &particles,
// 		     const float deltaT, int integrationMethod,
// 		     int plicCorrection, int vofCorrection,
// 		     int numSubSteps)
// {
//   vtkDataArray *velocityArray0 = velocityGrid[0]->GetCellData()->GetArray(0);
//   vtkDataArray *velocityArray1 = velocityGrid[1]->GetCellData()->GetArray(0);

//   int nodeRes[3];
//   velocityGrid[1]->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  
// #pragma omp parallel for  
//   for (int i = 0; i < particles.size(); ++i) {
//     particles[i] = rungeKutta4(particles[i], velocityGrid, velocityArray0, velocityArray1, cellRes,
// 			       0.0f, 1.0f, 0.0f, deltaT, numSubSteps);
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
		   std::vector<std::vector<int>> &neighbors,
		   const int procId)
{
  const int numDims = 3;
  const int numSides = 6;
  
  for (int i = 0; i < numDims; ++i) {

    if (localExtents[i*2+0] > globalExtents[i*2+0]) {
      for (int j = 0; j < allExtents.size()/numSides; ++j) {

	if (j == procId) continue;
	if (localExtents[i*2+0] <= allExtents[j*numSides+i*2+1] &&
	    localExtents[i*2+1] > allExtents[j*numSides+i*2+1] &&
	    localExtents[((i+1)%3)*2+0] <= allExtents[j*numSides+((i+1)%3)*2+1] &&
	    localExtents[((i+1)%3)*2+1] >= allExtents[j*numSides+((i+1)%3)*2+0] &&
	    localExtents[((i+2)%3)*2+0] <= allExtents[j*numSides+((i+2)%3)*2+1] &&
	    localExtents[((i+2)%3)*2+1] >= allExtents[j*numSides+((i+2)%3)*2+0]) {

	  neighbors[i*2+0].push_back(j);
	}
      }
    }
    if (localExtents[i*2+1] < globalExtents[i*2+1]) { 
      for (int j = 0; j < allExtents.size()/numSides; ++j) {

	if (j == procId) continue;
	if (localExtents[i*2+1] >= allExtents[j*numSides+i*2+0] &&
	    localExtents[i*2+0] < allExtents[j*numSides+i*2+0] &&
	    localExtents[((i+1)%3)*2+0] <= allExtents[j*numSides+((i+1)%3)*2+1] &&
	    localExtents[((i+1)%3)*2+1] >= allExtents[j*numSides+((i+1)%3)*2+0] &&
	    localExtents[((i+2)%3)*2+0] <= allExtents[j*numSides+((i+2)%3)*2+1] &&
	    localExtents[((i+2)%3)*2+1] >= allExtents[j*numSides+((i+2)%3)*2+0]) {

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
			 const int myExtent[6], const int globalExtent[6],
			 int cellRes[3], vtkFloatArray *labels,
			 std::vector<std::vector<float4> > &labelsToSend, int numGhosts)
{
  const int NUM_SIDES = 6;
  int slabs[NUM_SIDES][6] =
    {{0,numGhosts, 0,cellRes[1]-1, 0,cellRes[2]-1},
     {cellRes[0]-1-numGhosts,cellRes[0]-1, 0,cellRes[1]-1, 0,cellRes[2]-1},
     {0,cellRes[0]-1, 0,numGhosts, 0,cellRes[2]-1},
     {0,cellRes[0]-1, cellRes[1]-1-numGhosts,cellRes[1]-1, 0,cellRes[2]-1},
     {0,cellRes[0]-1, 0,cellRes[1]-1, 0,numGhosts},
     {0,cellRes[0]-1, 0,cellRes[1]-1, cellRes[2]-1-numGhosts,cellRes[2]-1}};

  if (myExtent[0] == globalExtent[0]) {
    slabs[0][0] = 0;
    slabs[0][1] = 1;
    slabs[2][0] = 0;
    slabs[3][0] = 0;
    slabs[4][0] = 0;
    slabs[5][0] = 0;
  }
  if (myExtent[1] == globalExtent[1]) {
    slabs[1][0] = cellRes[0]-2;
    slabs[1][1] = cellRes[0]-1;
    slabs[2][1] = cellRes[0]-1;
    slabs[3][1] = cellRes[0]-1;
    slabs[4][1] = cellRes[0]-1;
    slabs[5][1] = cellRes[0]-1;    
  }
  if (myExtent[2] == globalExtent[2]) {
    slabs[0][2] = 0;
    slabs[1][2] = 0;
    slabs[2][2] = 0;
    slabs[2][3] = 1;
    slabs[4][2] = 0;
    slabs[5][2] = 0;
  }
  if (myExtent[3] == globalExtent[3]) {
    slabs[0][3] = cellRes[1]-1;
    slabs[1][3] = cellRes[1]-1;
    slabs[3][2] = cellRes[1]-2;
    slabs[3][3] = cellRes[1]-1;
    slabs[4][3] = cellRes[1]-1;
    slabs[5][3] = cellRes[1]-1;    
  }
  if (myExtent[4] == globalExtent[4]) {
    slabs[0][4] = 0;
    slabs[1][4] = 0;
    slabs[2][4] = 0;
    slabs[3][4] = 0;
    slabs[4][4] = 0;
    slabs[4][5] = 1;
  }
  if (myExtent[5] == globalExtent[5]) {
    slabs[0][5] = cellRes[2]-1;
    slabs[1][5] = cellRes[2]-1;
    slabs[2][5] = cellRes[2]-1;
    slabs[3][5] = cellRes[2]-1;
    slabs[5][4] = cellRes[2]-2;
    slabs[5][5] = cellRes[2]-1;    
  }
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

	if (x >= myExtent[0] && x < myExtent[1] &&
	    y >= myExtent[2] && y < myExtent[3] &&
	    z >= myExtent[4] && z < myExtent[5]) {
	    
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

void calcLabelPoints(const std::vector<float> &labels,
		     const int labelsRange[2],
		     std::vector<std::vector<int>> &labelPoints)
{
  const int numPoints = labels.size();//->GetNumberOfTuples();
  const int numLabels = labelsRange[1] - labelsRange[0] + 1;
  for (int i = 0; i < numPoints; ++i) {

    int pointLabel = labels[i];//->GetComponent(i, 0);
    if (pointLabel < labelsRange[0] || pointLabel > labelsRange[1]) {
      continue;
    }    
    int pointLabelIdx = pointLabel - labelsRange[0];
    labelPoints[pointLabelIdx].push_back(i);
  }
}

void calcLabelExtents(const std::vector<float4> &points,
		      const std::vector<float> &labels,
		      const int labelsRange[2],
		      vtkRectilinearGrid *grid,
		      std::vector<std::array<int,6>> &labelBounds)
{
  const int numLabels = labelsRange[1] - labelsRange[0] + 1;
  std::vector<std::array<float,6>> labelBoundsTmp(numLabels);
  // initialize bounds
  for (int i = 0; i < numLabels; ++i) {
    labelBoundsTmp[i][0] = labelBoundsTmp[i][2] = labelBoundsTmp[i][4] = 
      std::numeric_limits<float>::max();
    labelBoundsTmp[i][1] = labelBoundsTmp[i][3] = labelBoundsTmp[i][5] = 
      -1.0f*std::numeric_limits<float>::max();
  }

  // compute physical bounds of each label
  const int numPoints = points.size();//->GetNumberOfPoints();
  for (int i = 0; i < numPoints; ++i) {
    
    int pointLabel = labels[i];//->GetComponent(i, 0);
    if (pointLabel < labelsRange[0] || pointLabel > labelsRange[1]) {
      continue;
    }    

    int pointLabelIdx = pointLabel - labelsRange[0];

    double p[3] = {points[i].x, points[i].y, points[i].z};
    // points->GetPoint(i, p);
    if (labelBoundsTmp[pointLabelIdx][0] > p[0]) labelBoundsTmp[pointLabelIdx][0] = p[0];
    if (labelBoundsTmp[pointLabelIdx][1] < p[0]) labelBoundsTmp[pointLabelIdx][1] = p[0];
    if (labelBoundsTmp[pointLabelIdx][2] > p[1]) labelBoundsTmp[pointLabelIdx][2] = p[1];
    if (labelBoundsTmp[pointLabelIdx][3] < p[1]) labelBoundsTmp[pointLabelIdx][3] = p[1];
    if (labelBoundsTmp[pointLabelIdx][4] > p[2]) labelBoundsTmp[pointLabelIdx][4] = p[2];
    if (labelBoundsTmp[pointLabelIdx][5] < p[2]) labelBoundsTmp[pointLabelIdx][5] = p[2];	
  }

  // transform physical bounds into grid bounds (index-based)
  for (int i = 0; i < numLabels; ++i) {
    
    double x0[3] = {labelBoundsTmp[i][0], labelBoundsTmp[i][2], labelBoundsTmp[i][4]};
    double x1[3] = {labelBoundsTmp[i][1], labelBoundsTmp[i][3], labelBoundsTmp[i][5]};
    int ijk0[3];
    int ijk1[3];
    double pcoords[3];

    grid->ComputeStructuredCoordinates(x0, ijk0, pcoords);
    grid->ComputeStructuredCoordinates(x1, ijk1, pcoords);

    labelBounds[i] = {ijk0[0],ijk1[0], ijk0[1],ijk1[1], ijk0[2],ijk1[2]};
  }
}

float4 computeNormal(const float *data, const int res[3],
		     const int idxCell[3], const double pcoords[3])
{
  int idxCell2[3] = {idxCell[0], idxCell[1], idxCell[2]};

  idxCell2[0] -= 1;
  float left = interpolateSca(data, res,  idxCell2, pcoords);
  idxCell2[0] = idxCell[0];

  idxCell2[0] += 1;
  float right = interpolateSca(data, res,  idxCell2, pcoords);
  idxCell2[0] = idxCell[0];

  idxCell2[1] -= 1;
  float bottom = interpolateSca(data, res,  idxCell2, pcoords);
  idxCell2[1] = idxCell[1];

  idxCell2[1] += 1;
  float top = interpolateSca(data, res,  idxCell2, pcoords);
  idxCell2[1] = idxCell[1];

  idxCell2[2] -= 1;
  float back = interpolateSca(data, res,  idxCell2, pcoords);
  idxCell2[2] = idxCell[2];

  idxCell2[2] += 1;
  float front = interpolateSca(data, res,  idxCell2, pcoords);
  idxCell2[2] = idxCell[2];

  float3 diff = make_float3(right-left, top-bottom, front-back);
  if (length(diff) > 0.0f) {
    diff = normalize(diff);
  }
  diff = diff*-1.0f;
  
  return make_float4(diff.x, diff.y, diff.z, 0.0f);
}

void computeNormals(vtkRectilinearGrid *grid,
		    const float *data,
		    std::vector<float4>::iterator itv,
		    std::vector<float4>::iterator itve,
		    std::vector<float4> &normals)
{
  int index;
  int nodeRes[3];
  grid->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  
  for (; itv != itve; ++itv) {
    double p[3] = {itv->x, itv->y, itv->z};
    int ijk[3];
    double pcoords[3];

    grid->ComputeStructuredCoordinates(p, ijk, pcoords);
    float4 n = computeNormal(data, nodeRes, ijk, pcoords);
    normals.push_back(n);
  }
}


void computeNormals(std::vector<float4> &vertices,
		    std::vector<int> &indices,
		    int begin, int end,
		    std::vector<float4> &normals)
{
  std::map<int,std::vector<int>> neighborVertices;
  neighborVertices.clear();

  for (int i = begin; i < end; i+=3) {
    for (int j = 0; j < 3; ++j) {
      int id0 = indices[i*3+j];
      int id1 = indices[i*3+((j+1)%3)];
      int id2 = indices[i*3+((j+2)%3)];
      if (neighborVertices.find(id0) == neighborVertices.end()) {
	std::vector<int> ve(2);
	ve[0] = id1;
	ve[1] = id2;
	neighborVertices[id0] = ve;
      }
      else {
	neighborVertices[id0].push_back(id1);
	neighborVertices[id0].push_back(id2);
      }
    }
  }

  // int numNormals = 0;
  int numNormals = normals.size();

  std::map<int,std::vector<int>>::iterator it = neighborVertices.begin();
  for (; it != neighborVertices.end(); ++it) {
    std::vector<int> &neighbors = it->second;
    float3 normal = make_float3(0.0f);
    float3 v0 = make_float3(vertices[it->first]);
    for (int i = 0; i < neighbors.size(); i+=2) {

      float3 v1 = make_float3(vertices[neighbors[i]]);
      float3 v2 = make_float3(vertices[neighbors[i+1]]);
      float3 e1 = v1 - v0;
      float3 e2 = v2 - v0;
      normal += cross(e1,e2);
    }
    normals.push_back(make_float4(normalize(normal)));
    // ++numNormals;
  }

  // smooth normals
  // std::vector<float4> normalsTmp(numNormals);
  // for (int i = 0; i < numNormals; ++i)
  //   normalsTmp[i] = make_float4(0.0f);
  std::vector<float4> normalsTmp(0);
  
  for (; it != neighborVertices.end(); ++it) {
    std::vector<int> &neighbors = it->second;
    float3 n0 = make_float3(normals[it->first]);
    float3 v0 = make_float3(vertices[it->first]);

    float3 normal = make_float3(0.0f);
    for (int i = 0; i < neighbors.size(); i+=2) {

      float3 v1 = make_float3(vertices[neighbors[i]]);
      float3 v2 = make_float3(vertices[neighbors[i+1]]);
      float3 n1 = make_float3(normals[neighbors[i]]);
      float3 n2 = make_float3(normals[neighbors[i+1]]);
      float3 e1 = v1 - v0;
      float3 e2 = v2 - v0;
      float area = length(cross(e1,e2));
      float3 sumn = (n0+n1+n2)*area;
      normal += sumn;
    }
    normalsTmp.push_back(make_float4(normalize(normal)));
  }

  for (int i = 0; i < normalsTmp.size(); ++i) {
    normals[numNormals+i] = normalsTmp[i];
  }
}

void resamplePointsOnGrid(const std::vector<int> &labelPoints,
			  const std::vector<float4> &points,
			  vtkRectilinearGrid *subGrid,
			  const int subNodeRes[3],
			  std::vector<float> &field)
{
  for (int j = 0; j < labelPoints.size(); ++j) {

    float4 px = points[labelPoints[j]];
    double x[3] = {px.x, px.y, px.z};
    int ijk[3];
    double pcoords[3];    
    int inside = subGrid->ComputeStructuredCoordinates(x, ijk, pcoords);

    // pcoords[0] += 0.5;
    // pcoords[1] += 0.5;
    // pcoords[2] += 0.5;

    if (!inside) continue;
    
    int ids[8] =
      {ijk[0]   +  ijk[1]*subNodeRes[0]    +  ijk[2]*subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 +  ijk[1]*subNodeRes[0]    +  ijk[2]*subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 + (ijk[1]+1)*subNodeRes[0] +  ijk[2]*subNodeRes[0]*subNodeRes[1],
       ijk[0]   + (ijk[1]+1)*subNodeRes[0] +  ijk[2]*subNodeRes[0]*subNodeRes[1],
       ijk[0]   +  ijk[1]*subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 +  ijk[1]*subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 + (ijk[1]+1)*subNodeRes[0] + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
       ijk[0]   + (ijk[1]+1)*subNodeRes[0] + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1]};

    field[ids[0]] += (1.0f-pcoords[0])*(1.0f-pcoords[1])*(1.0f-pcoords[2]);
    field[ids[1]] += (     pcoords[0])*(1.0f-pcoords[1])*(1.0f-pcoords[2]);
    field[ids[2]] += (     pcoords[0])*(     pcoords[1])*(1.0f-pcoords[2]);
    field[ids[3]] += (1.0f-pcoords[0])*(     pcoords[1])*(1.0f-pcoords[2]);
    field[ids[4]] += (1.0f-pcoords[0])*(1.0f-pcoords[1])*(     pcoords[2]);
    field[ids[5]] += (     pcoords[0])*(1.0f-pcoords[1])*(     pcoords[2]);
    field[ids[6]] += (     pcoords[0])*(     pcoords[1])*(     pcoords[2]);
    field[ids[7]] += (1.0f-pcoords[0])*(     pcoords[1])*(     pcoords[2]);
  }
}

void resamplePointsOnGrid2(const std::vector<int> &labelPoints,
			  const std::vector<float> &labels,
			  const std::vector<float4> &points,
			  vtkRectilinearGrid *subGrid,
			  const int subNodeRes[3],
			  std::vector<float> &field)
{
  for (int j = 0; j < labelPoints.size(); ++j) {

    float4 px = points[labelPoints[j]];
    float label = labels[labelPoints[j]];
    double x[3] = {px.x, px.y, px.z};
    int ijk[3];
    double pcoords[3];    
    int inside = subGrid->ComputeStructuredCoordinates(x, ijk, pcoords);

    if (!inside) continue;
    
    int ids[8] =
      {ijk[0]   +  ijk[1]   *subNodeRes[0]    +  ijk[2]   *subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 +  ijk[1]   *subNodeRes[0]    +  ijk[2]   *subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 + (ijk[1]+1)*subNodeRes[0]    +  ijk[2]   *subNodeRes[0]*subNodeRes[1],
       ijk[0]   + (ijk[1]+1)*subNodeRes[0]    +  ijk[2]   *subNodeRes[0]*subNodeRes[1],
       ijk[0]   +  ijk[1]   *subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 +  ijk[1]   *subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
       ijk[0]+1 + (ijk[1]+1)*subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
       ijk[0]   + (ijk[1]+1)*subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1]};

    field[ids[0]] += (label+1)*(1.0f-pcoords[0])*(1.0f-pcoords[1])*(1.0f-pcoords[2]);
    field[ids[1]] += (label+1)*(     pcoords[0])*(1.0f-pcoords[1])*(1.0f-pcoords[2]);
    field[ids[2]] += (label+1)*(     pcoords[0])*(     pcoords[1])*(1.0f-pcoords[2]);
    field[ids[3]] += (label+1)*(1.0f-pcoords[0])*(     pcoords[1])*(1.0f-pcoords[2]);
    field[ids[4]] += (label+1)*(1.0f-pcoords[0])*(1.0f-pcoords[1])*(     pcoords[2]);
    field[ids[5]] += (label+1)*(     pcoords[0])*(1.0f-pcoords[1])*(     pcoords[2]);
    field[ids[6]] += (label+1)*(     pcoords[0])*(     pcoords[1])*(     pcoords[2]);
    field[ids[7]] += (label+1)*(1.0f-pcoords[0])*(     pcoords[1])*(     pcoords[2]);
  }
}


void generateCoords(vtkDataArray *coords,
		    const int subNodeRes,
		    const int subone,
		    const int r,
		    const int ijk0, const int ijk1,
		    vtkFloatArray **subcoords)
{
  (*subcoords) = vtkFloatArray::New();
  (*subcoords)->SetNumberOfComponents(1);
  (*subcoords)->SetNumberOfTuples(subNodeRes);

  float xprev = coords->GetComponent(ijk0, 0);    
  int ires = (subNodeRes + (r-1))/r;   // subNodeRes = subNodeRes*r - (r-1);

  for (int j = 0; j < ires-1; ++j) {
    float x = coords->GetComponent(ijk0+j+1, 0);
    float dx = (x - xprev)/r;
	
    for (int k = 0; k < r; ++k) {
      (*subcoords)->SetValue(j*r+k, xprev+k*dx);
    }
    xprev = x;
  }
  (*subcoords)->SetValue(subNodeRes-1, coords->GetComponent(ijk1+1,0)); 
}

void generateBoundary(const std::vector<float4> &points,
		      const std::vector<float> &labels,
		      vtkRectilinearGrid *grid,
		      const int refinement,
		      const int localExtentNoGhosts[6],
		      const int localExtent[6],
		      int &vertexID,
		      std::vector<int> &labelOffsets,
		      std::vector<float4> &vertices,
		      std::vector<float4> &normals,
		      std::vector<int> &indices)
{
  const float isoValue = 0.75f;
  vtkDataArray *coords[3] = {grid->GetXCoordinates(), 
			     grid->GetYCoordinates(), 
			     grid->GetZCoordinates()};
  int nodeRes[3];
  grid->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

  const int numPoints = points.size();
  const int boundarySize = 2; // for correct normal computation
  double range[2] = {-1, *std::max_element(labels.begin(), labels.end())};

  const int numLabels = std::ceil(range[1] - range[0] + 1.0f);
  const int labelsRange[2] = {range[0], range[1]};

  // stores indices of all points with a given label
  std::vector<std::vector<int>> labelPoints(numLabels);
  for (int i = 0; i < numLabels; ++i) {
    labelPoints[i].clear();
  }
  calcLabelPoints(labels, labelsRange, labelPoints);

  // stores the index bounds of particles with a given label
  // the bounds are given in coordinates relative to the local extent
  std::vector<std::array<int,6>> labelExtents(numLabels);
  calcLabelExtents(points, labels, labelsRange, grid, labelExtents);

  labelOffsets.resize(numLabels+1,0);
  
  for (int i = 0; i < numLabels; ++i) {

    if (labelPoints[i].size() == 0) {
      labelOffsets[i+1] = vertices.size();
      continue;
    }

    const std::array<int,6> &labelExtent = labelExtents[i];
    const std::array<int,6> nodeExtent = {labelExtent[0]-1,labelExtent[1]+1,
    					  labelExtent[2]-1,labelExtent[3]+1,
    					  labelExtent[4]-1,labelExtent[5]+1};
    int subNodeRes[3] = {nodeExtent[1]-nodeExtent[0]+1,
			 nodeExtent[3]-nodeExtent[2]+1,
			 nodeExtent[5]-nodeExtent[4]+1};
    
    vtkFloatArray *subcoords[3];

    for (int j = 0; j < 3; ++j) {
      subcoords[j] = vtkFloatArray::New();
      subcoords[j]->SetNumberOfComponents(1);
      subcoords[j]->SetNumberOfTuples(subNodeRes[j]);

      int beg = 0;
      if (nodeExtent[j*2] < 0) {
	float x0 = coords[j]->GetComponent(0,0);
	float x1 = coords[j]->GetComponent(1,0);
	subcoords[j]->SetComponent(0,0,x0-(x1-x0)/2.0f);
	beg = 1;
      }
      int end = subNodeRes[j];
      if (nodeExtent[j*2] + end > nodeRes[j]-1) {
	float x0 = coords[j]->GetComponent(nodeRes[j]-2,0);
	float x1 = coords[j]->GetComponent(nodeRes[j]-1,0);
	subcoords[j]->SetComponent(subNodeRes[j]-1,0,x1+(x1-x0)/2.0f);
	end = subNodeRes[j]-1;	
      }

      for (int c = beg; c < end; ++c) {
	float x0 = coords[j]->GetComponent(nodeExtent[j*2]+c,0);
	float x1 = coords[j]->GetComponent(nodeExtent[j*2]+c+1,0);
	subcoords[j]->SetComponent(c,0,(x1+x0)/2.0f);
      }
    }

    const int r = std::pow(2,refinement);
    if (refinement > 0) {

      for (int j = 0; j < 3; ++j) {

	int prev_snr = subNodeRes[j];
	subNodeRes[j] = (subNodeRes[j]-1)*r;

	vtkFloatArray *subcoord_tmp = vtkFloatArray::New();
	subcoord_tmp->SetNumberOfComponents(1);
	subcoord_tmp->SetNumberOfTuples(subNodeRes[j]);

	for (int c = 0; c < prev_snr-1; ++c) {
	  float dx = subcoords[j]->GetComponent(c+1,0)-subcoords[j]->GetComponent(c,0);
	  dx /= r;
	  float org = subcoords[j]->GetComponent(c,0) + dx/2.0f;
	  for (int s = 0; s < r; ++s) {
	    subcoord_tmp->SetComponent(c*r+s,0,org+dx*s);
	  }	  
	}
      	subcoords[j]->Delete();
      	subcoords[j] = subcoord_tmp;
      }
    }

    vtkRectilinearGrid *subGrid = vtkRectilinearGrid::New();    
    subGrid->SetDimensions(subcoords[0]->GetNumberOfTuples(),
			   subcoords[1]->GetNumberOfTuples(),
			   subcoords[2]->GetNumberOfTuples());   
    subGrid->SetXCoordinates(subcoords[0]);
    subGrid->SetYCoordinates(subcoords[1]);
    subGrid->SetZCoordinates(subcoords[2]);

    const int numElements = subNodeRes[0]*subNodeRes[1]*subNodeRes[2];
    std::vector<float> field(numElements, 0.0f);

    resamplePointsOnGrid(labelPoints[i], points, subGrid, subNodeRes, field);    

    int subGridExtent[6] = {0, subNodeRes[0]-1,
			    0, subNodeRes[1]-1,
			    0, subNodeRes[2]-1};

    int df;
    df = labelExtent[0]+localExtent[0] - localExtentNoGhosts[0];
    if (df < 0) subGridExtent[0] += r*std::min(-df,boundarySize);

    df = labelExtent[1]+localExtent[0] - localExtentNoGhosts[1]; 
    if (df > 0) subGridExtent[1] -= r*std::min(df,boundarySize);

    df = labelExtent[2]+localExtent[2] - localExtentNoGhosts[2];
    if (df < 0) subGridExtent[2] += r*std::min(-df,boundarySize);

    df = labelExtent[3]+localExtent[2] - localExtentNoGhosts[3];
    if (df > 0) subGridExtent[3] -= r*std::min(df,boundarySize);

    df = labelExtent[4]+localExtent[4] - localExtentNoGhosts[4];
    if (df < 0) subGridExtent[4] += r*std::min(-df,boundarySize);

    df = labelExtent[5]+localExtent[4] - localExtentNoGhosts[5]; 
    if (df > 0) subGridExtent[5] -= r*std::min(df,boundarySize);

    int numVertsPrev = vertices.size();
    int numIndicesPrev = indices.size();
    extractSurface(field.data(), subNodeRes, subcoords, subGridExtent,
		   isoValue, indices, vertices, vertexID);

    computeNormals(subGrid, field.data(), vertices.begin()+numVertsPrev, vertices.end(), normals);
    
    subGrid->Delete();

    labelOffsets[i+1] = vertices.size();
  }
}

void calcLabelExtents(const std::vector<float4> &points,
		      const std::vector<int> &labelPoints,
		      vtkRectilinearGrid *grid,
		      std::array<int,6> &labelBounds)
{
  std::array<float,6> labelBoundsTmp;

  // initialize bounds
  labelBoundsTmp[0] = labelBoundsTmp[2] = labelBoundsTmp[4] = std::numeric_limits<float>::max();
  labelBoundsTmp[1] = labelBoundsTmp[3] = labelBoundsTmp[5] = -1.0f*std::numeric_limits<float>::max();

  // compute physical bounds of each label
  const int numPoints = labelPoints.size();
  for (int j = 0; j < numPoints; ++j) {

    int i = labelPoints[j];
    double p[3] = {points[i].x, points[i].y, points[i].z};
    // points->GetPoint(i, p);
    if (labelBoundsTmp[0] > p[0]) labelBoundsTmp[0] = p[0];
    if (labelBoundsTmp[1] < p[0]) labelBoundsTmp[1] = p[0];
    if (labelBoundsTmp[2] > p[1]) labelBoundsTmp[2] = p[1];
    if (labelBoundsTmp[3] < p[1]) labelBoundsTmp[3] = p[1];
    if (labelBoundsTmp[4] > p[2]) labelBoundsTmp[4] = p[2];
    if (labelBoundsTmp[5] < p[2]) labelBoundsTmp[5] = p[2];	
  }

  double x0[3] = {labelBoundsTmp[0], labelBoundsTmp[2], labelBoundsTmp[4]};
  double x1[3] = {labelBoundsTmp[1], labelBoundsTmp[3], labelBoundsTmp[5]};
  int ijk0[3];
  int ijk1[3];
  double pcoords[3];
  
  grid->ComputeStructuredCoordinates(x0, ijk0, pcoords);
  grid->ComputeStructuredCoordinates(x1, ijk1, pcoords);

  labelBounds = {ijk0[0],ijk1[0], ijk0[1],ijk1[1], ijk0[2],ijk1[2]};
}

void generatePairs(std::unordered_set<float> labels,
		   std::vector<std::pair<float,float>> &pairs)
{
  std::vector<float> labelsVec;
  labelsVec.clear();
  for (auto it = labels.begin(); it != labels.end(); ++it) {
    labelsVec.push_back(*it);
  }

  std::sort(labelsVec.begin(), labelsVec.end());
  
  for (int i = 0; i < labelsVec.size(); ++i) {
    for (int j = i+1; j < labelsVec.size(); ++j) {
      pairs.push_back(std::pair<float,float>(labelsVec[i],labelsVec[j]));
    }
  }
}

void generateBoundary(const std::vector<float4> &points,
		      const std::vector<float> &labels,		      
		      vtkRectilinearGrid *grid,
		      const int refinement,
		      const int localExtentNoGhosts[6],
		      const int localExtent[6],
		      int &vertexID,
		      std::vector<float4> &vertices,
		      std::vector<float4> &normals,
		      std::vector<int> &indices,
		      std::vector<std::vector<int>> &prevLabelPoints)
{
  {
    long long numPoints = points.size();
    long long numPrevPoints = 0;
    for (const auto  &it : prevLabelPoints) {
      numPrevPoints += it.size();
    }
    if (numPrevPoints != 0 && numPoints != numPrevPoints) {
      std::cout << "ERROR: the number of points doesn't match: "
		<< numPoints << " - " << numPrevPoints << std::endl;
    }
  }
  
  if (points.empty()) {
    return;
  }
  vtkDataArray *coords[3] = {grid->GetXCoordinates(), 
			     grid->GetYCoordinates(), 
			     grid->GetZCoordinates()};
  int nodeRes[3];
  grid->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

  const int numPoints = points.size();
  const int boundarySize = 2; // for correct normal computation
  double range[2] = {-1, *std::max_element(labels.begin(), labels.end())};

  const int numLabels = std::ceil(range[1] - range[0] + 1.0f);
  const int labelsRange[2] = {range[0], range[1]};

  // stores indices of all points with a given label
  std::vector<std::vector<int>> labelPoints(numLabels);
  for (int i = 0; i < numLabels; ++i) {
    labelPoints[i].clear();
  }
  calcLabelPoints(labels, labelsRange, labelPoints);

  for (int i = 0; i < prevLabelPoints.size(); ++i) {

    //---------------------------------------------------
    std::unordered_set<float> unique_labels;
    unique_labels.clear();
    for (int j = 0; j < prevLabelPoints[i].size(); ++j) {
      // if (labels[prevLabelPoints[i][j]] != -1) {
	unique_labels.insert(labels[prevLabelPoints[i][j]]);
      // }
    }
    if (unique_labels.size() < 2) {
      continue;
    }
    //---------------------------------------------------
    
    std::array<int,6> labelExtent;
    calcLabelExtents(points, prevLabelPoints[i], grid, labelExtent);

    const std::array<int,6> nodeExtent = {labelExtent[0]-1,labelExtent[1]+1,
    					  labelExtent[2]-1,labelExtent[3]+1,
    					  labelExtent[4]-1,labelExtent[5]+1};
    int subNodeRes[3] = {nodeExtent[1]-nodeExtent[0]+1,
			 nodeExtent[3]-nodeExtent[2]+1,
			 nodeExtent[5]-nodeExtent[4]+1};
    
    vtkFloatArray *subcoords[3];

    for (int j = 0; j < 3; ++j) {
      subcoords[j] = vtkFloatArray::New();
      subcoords[j]->SetNumberOfComponents(1);
      subcoords[j]->SetNumberOfTuples(subNodeRes[j]);

      int beg = 0;
      if (nodeExtent[j*2] < 0) {
	float x0 = coords[j]->GetComponent(0,0);
	float x1 = coords[j]->GetComponent(1,0);
	subcoords[j]->SetComponent(0,0,x0-(x1-x0)/2.0f);
	beg = 1;
      }
      int end = subNodeRes[j];
      if (nodeExtent[j*2] + end > nodeRes[j]-1) {
	float x0 = coords[j]->GetComponent(nodeRes[j]-2,0);
	float x1 = coords[j]->GetComponent(nodeRes[j]-1,0);
	subcoords[j]->SetComponent(subNodeRes[j]-1,0,x1+(x1-x0)/2.0f);
	end = subNodeRes[j]-1;	
      }

      for (int c = beg; c < end; ++c) {
	float x0 = coords[j]->GetComponent(nodeExtent[j*2]+c,0);
	float x1 = coords[j]->GetComponent(nodeExtent[j*2]+c+1,0);
	subcoords[j]->SetComponent(c,0,(x1+x0)/2.0f);
      }
    }

    const int r = std::pow(2,refinement);
    if (refinement > 0) {

      for (int j = 0; j < 3; ++j) {

	int prev_snr = subNodeRes[j];
	subNodeRes[j] = (subNodeRes[j]-1)*r;

	vtkFloatArray *subcoord_tmp = vtkFloatArray::New();
	subcoord_tmp->SetNumberOfComponents(1);
	subcoord_tmp->SetNumberOfTuples(subNodeRes[j]);

	for (int c = 0; c < prev_snr-1; ++c) {
	  float dx = subcoords[j]->GetComponent(c+1,0)-subcoords[j]->GetComponent(c,0);
	  dx /= r;
	  float org = subcoords[j]->GetComponent(c,0) + dx/2.0f;
	  for (int s = 0; s < r; ++s) {
	    subcoord_tmp->SetComponent(c*r+s,0,org+dx*s);
	  }	  
	}
      	subcoords[j]->Delete();
      	subcoords[j] = subcoord_tmp;
      }
    }

    vtkRectilinearGrid *subGrid = vtkRectilinearGrid::New();    
    subGrid->SetDimensions(subcoords[0]->GetNumberOfTuples(),
			   subcoords[1]->GetNumberOfTuples(),
			   subcoords[2]->GetNumberOfTuples());   
    subGrid->SetXCoordinates(subcoords[0]);
    subGrid->SetYCoordinates(subcoords[1]);
    subGrid->SetZCoordinates(subcoords[2]);

    const int numElements = subNodeRes[0]*subNodeRes[1]*subNodeRes[2];
    std::vector<float> field(numElements, -1.0f);

    resamplePointsOnGrid2(prevLabelPoints[i], labels, points, subGrid, subNodeRes, field);    

    int subGridExtent[6] = {0, subNodeRes[0]-1,
			    0, subNodeRes[1]-1,
			    0, subNodeRes[2]-1};

    int df;
    df = labelExtent[0]+localExtent[0] - localExtentNoGhosts[0];
    if (df < 0) subGridExtent[0] += r*std::min(-df,boundarySize);

    df = labelExtent[1]+localExtent[0] - localExtentNoGhosts[1]; 
    if (df > 0) subGridExtent[1] -= r*std::min(df,boundarySize);

    df = labelExtent[2]+localExtent[2] - localExtentNoGhosts[2];
    if (df < 0) subGridExtent[2] += r*std::min(-df,boundarySize);

    df = labelExtent[3]+localExtent[2] - localExtentNoGhosts[3];
    if (df > 0) subGridExtent[3] -= r*std::min(df,boundarySize);

    df = labelExtent[4]+localExtent[4] - localExtentNoGhosts[4];
    if (df < 0) subGridExtent[4] += r*std::min(-df,boundarySize);

    df = labelExtent[5]+localExtent[4] - localExtentNoGhosts[5]; 
    if (df > 0) subGridExtent[5] -= r*std::min(df,boundarySize);

    const float isoValue = 0.75f;
    int numVertsPrev = vertices.size();
    int numIndicesPrev = indices.size();


    std::vector<std::pair<float,float>> combs;
    combs.clear();

    generatePairs(unique_labels, combs);

    for (int j = 0; j < combs.size(); ++j) {
    
      extractSurface2(field.data(), subNodeRes, subcoords, subGridExtent,
		      isoValue, indices, vertices, vertexID, combs[j].first, combs[j].second);
    }

    computeNormals(subGrid, field.data(), vertices.begin()+numVertsPrev, vertices.end(), normals);
    
    subGrid->Delete();
  }

  prevLabelPoints = labelPoints;
}

// void generateBoundary(const std::vector<float4> &points,
// 		      const std::vector<float> &labels,		      
// 		      vtkRectilinearGrid *grid,
// 		      const int refinement,
// 		      const int localExtentNoGhosts[6],
// 		      const int localExtent[6],
// 		      int &vertexID,
// 		      std::vector<float4> &vertices,
// 		      std::vector<float4> &normals,
// 		      std::vector<int> &indices,
// 		      std::vector<std::vector<int>> &prevLabelPoints)
// {
//   const float isoValue = 0.75f;
//   vtkDataArray *coords[3] = {grid->GetXCoordinates(), 
// 			     grid->GetYCoordinates(), 
// 			     grid->GetZCoordinates()};
//   int nodeRes[3];
//   grid->GetDimensions(nodeRes);
//   int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};

//   const int numPoints = points.size();
//   const int boundarySize = 2; // for correct normal computation
//   double range[2] = {-1, *std::max_element(labels.begin(), labels.end())};

//   const int numLabels = std::ceil(range[1] - range[0] + 1.0f);
//   const int labelsRange[2] = {range[0], range[1]};

//   // stores indices of all points with a given label
//   std::vector<std::vector<int>> labelPoints(numLabels);
//   for (int i = 0; i < numLabels; ++i) {
//     labelPoints[i].clear();
//   }
//   calcLabelPoints(labels, labelsRange, labelPoints);


//   /*
//     prev    |-------|------|
//     curr    |--|------|----|
//    */
//   // // new ---------------------------------------------------------------------
//   // std::vector<bool> updateLabels(numLabels, false);

//   // if (!prevLabelPoints.empty()) {

//   //   int prev_beg = 0;
//   //   int prev_end = prevLabelPoints[0].size();
//   //   int curr_beg = 0;
//   //   int curr_end = labelPoints[0].size();
//   //   int curr_i = 0;
//   //   int prev_i = 0;

//   //   while (curr_i < numLabels) {
//   //     if (curr_end < prev_end) {
//   // 	while (curr_end < prev_end) {
	   
//   // 	  updateLabels[curr_i] = true;

//   // 	  curr_i++; // 1
//   // 	  curr_end += labelPoints[curr_i].size();
//   // 	}
//   // 	updateLabels[curr_i] = true; // 1
//   //     }
//   //     if (curr_end == prev_end) {

//   // 	prev_beg = prev_end;
//   // 	curr_beg = curr_end;	  

//   // 	prev_i++;
//   // 	prev_end += prevLabelPoints[prev_i].size();
//   // 	curr_i++;
//   // 	curr_end += labelPoints[curr_i].size();
//   //     }
//   //     if (curr_end > prev_end) {
//   // 	while (curr_end > prev_end) {
//   // 	  prev_i++; // 1
//   // 	  prev_end += prevLabelPoints[prev_i].size(); // 1
//   // 	}
//   //     }
//   //   }
//   // }
//   //~new ---------------------------------------------------------------------
  
//   // stores the index bounds of particles with a given label
//   // the bounds are given in coordinates relative to the local extent
//   std::vector<std::array<int,6>> labelExtents(numLabels);
//   calcLabelExtents(points, labels, labelsRange, grid, labelExtents);

//   for (int i = 0; i < numLabels; ++i) {

//     // if (!updateLabels[i]) {
//     //   continue;
//     // }

//     const std::array<int,6> &labelExtent = labelExtents[i];
//     const std::array<int,6> nodeExtent = {labelExtent[0]-1,labelExtent[1]+1,
//     					  labelExtent[2]-1,labelExtent[3]+1,
//     					  labelExtent[4]-1,labelExtent[5]+1};
//     int subNodeRes[3] = {nodeExtent[1]-nodeExtent[0]+1,
// 			 nodeExtent[3]-nodeExtent[2]+1,
// 			 nodeExtent[5]-nodeExtent[4]+1};
    
//     vtkFloatArray *subcoords[3];

//     for (int j = 0; j < 3; ++j) {
//       subcoords[j] = vtkFloatArray::New();
//       subcoords[j]->SetNumberOfComponents(1);
//       subcoords[j]->SetNumberOfTuples(subNodeRes[j]);

//       int beg = 0;
//       if (nodeExtent[j*2] < 0) {
// 	float x0 = coords[j]->GetComponent(0,0);
// 	float x1 = coords[j]->GetComponent(1,0);
// 	subcoords[j]->SetComponent(0,0,x0-(x1-x0)/2.0f);
// 	beg = 1;
//       }
//       int end = subNodeRes[j];
//       if (nodeExtent[j*2] + end > nodeRes[j]-1) {
// 	float x0 = coords[j]->GetComponent(nodeRes[j]-2,0);
// 	float x1 = coords[j]->GetComponent(nodeRes[j]-1,0);
// 	subcoords[j]->SetComponent(subNodeRes[j]-1,0,x1+(x1-x0)/2.0f);
// 	end = subNodeRes[j]-1;	
//       }

//       for (int c = beg; c < end; ++c) {
// 	float x0 = coords[j]->GetComponent(nodeExtent[j*2]+c,0);
// 	float x1 = coords[j]->GetComponent(nodeExtent[j*2]+c+1,0);
// 	subcoords[j]->SetComponent(c,0,(x1+x0)/2.0f);
//       }
//     }

//     const int r = std::pow(2,refinement);
//     if (refinement > 0) {

//       for (int j = 0; j < 3; ++j) {

// 	int prev_snr = subNodeRes[j];
// 	subNodeRes[j] = (subNodeRes[j]-1)*r;

// 	vtkFloatArray *subcoord_tmp = vtkFloatArray::New();
// 	subcoord_tmp->SetNumberOfComponents(1);
// 	subcoord_tmp->SetNumberOfTuples(subNodeRes[j]);

// 	for (int c = 0; c < prev_snr-1; ++c) {
// 	  float dx = subcoords[j]->GetComponent(c+1,0)-subcoords[j]->GetComponent(c,0);
// 	  dx /= r;
// 	  float org = subcoords[j]->GetComponent(c,0) + dx/2.0f;
// 	  for (int s = 0; s < r; ++s) {
// 	    subcoord_tmp->SetComponent(c*r+s,0,org+dx*s);
// 	  }	  
// 	}
//       	subcoords[j]->Delete();
//       	subcoords[j] = subcoord_tmp;
//       }
//     }

//     vtkRectilinearGrid *subGrid = vtkRectilinearGrid::New();    
//     subGrid->SetDimensions(subcoords[0]->GetNumberOfTuples(),
// 			   subcoords[1]->GetNumberOfTuples(),
// 			   subcoords[2]->GetNumberOfTuples());   
//     subGrid->SetXCoordinates(subcoords[0]);
//     subGrid->SetYCoordinates(subcoords[1]);
//     subGrid->SetZCoordinates(subcoords[2]);

//     const int numElements = subNodeRes[0]*subNodeRes[1]*subNodeRes[2];
//     std::vector<float> field(numElements, 0.0f);

//     resamplePointsOnGrid(labelPoints[i], points, subGrid, subNodeRes, field);    

//     int subGridExtent[6] = {0, subNodeRes[0]-1,
// 			    0, subNodeRes[1]-1,
// 			    0, subNodeRes[2]-1};

//     int df;
//     df = labelExtent[0]+localExtent[0] - localExtentNoGhosts[0];
//     if (df < 0) subGridExtent[0] += r*std::min(-df,boundarySize);

//     df = labelExtent[1]+localExtent[0] - localExtentNoGhosts[1]; 
//     if (df > 0) subGridExtent[1] -= r*std::min(df,boundarySize);

//     df = labelExtent[2]+localExtent[2] - localExtentNoGhosts[2];
//     if (df < 0) subGridExtent[2] += r*std::min(-df,boundarySize);

//     df = labelExtent[3]+localExtent[2] - localExtentNoGhosts[3];
//     if (df > 0) subGridExtent[3] -= r*std::min(df,boundarySize);

//     df = labelExtent[4]+localExtent[4] - localExtentNoGhosts[4];
//     if (df < 0) subGridExtent[4] += r*std::min(-df,boundarySize);

//     df = labelExtent[5]+localExtent[4] - localExtentNoGhosts[5]; 
//     if (df > 0) subGridExtent[5] -= r*std::min(df,boundarySize);

//     int numVertsPrev = vertices.size();
//     int numIndicesPrev = indices.size();
//     extractSurface(field.data(), subNodeRes, subcoords, subGridExtent,
// 		   isoValue, indices, vertices, vertexID);

//     computeNormals(subGrid, field.data(), vertices.begin()+numVertsPrev, vertices.end(), normals);
    
//     subGrid->Delete();
//   }

//   prevLabelPoints = labelPoints;
// }

void writeData(vtkPolyData *data, const int blockId,
	       const int processId, const std::string path)
{
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(data);
  std::string outFile = path;
  outFile += std::to_string(blockId) + std::string("_") + std::to_string(processId);
  outFile += ".vtp";
  writer->SetFileName(outFile.c_str());
  writer->SetDataModeToBinary();
  writer->Write();
}
void writeData(vtkRectilinearGrid *data, const int blockId,
	       const int processId, const std::string path)
{
  vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
  writer->SetInputData(data);
  std::string outFile = path;
  outFile += std::to_string(blockId) + std::string("_") + std::to_string(processId);
  outFile += ".vtr";
  writer->SetFileName(outFile.c_str());
  writer->SetDataModeToBinary();
  writer->Write();
}

