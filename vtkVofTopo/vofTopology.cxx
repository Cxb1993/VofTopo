#include "vofTopology.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkShortArray.h"
#include "vtkCellArray.h"
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <set>
#include <cmath>
#include <array>

#include "marchingCubes_cpu.h"

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


void initVelocities(vtkRectilinearGrid *velocity,
		    std::vector<float4> &particles,
		    std::vector<float4> &velocities)
{
  int nodeRes[3];
  velocity->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  double x[3];
  int ijk[3];
  double pcoords[3];
  int index;
  vtkDataArray *velocityArray = 
    velocity->GetCellData()->GetArray("Data", index);
      // velocity->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
  if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }

  std::vector<float4>::iterator itp = particles.begin();
  std::vector<float4>::iterator itv = velocities.begin();
  for (; itp != particles.end() && itv != velocities.end(); ++itp, ++itv) {

    x[0] = itp->x;
    x[1] = itp->y;
    x[2] = itp->z;        
    velocity->ComputeStructuredCoordinates(x, ijk, pcoords);
    *itv = make_float4(interpolateVec(velocityArray, cellRes, ijk, pcoords),0.0f);
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

void generateSeedPoints(vtkRectilinearGrid *velocityGrid,
			int refinement,
			vtkPoints *points,
			int globalExtent[6],
			int numGhostLevels)
{
  vtkDataArray *coords[3] = {velocityGrid->GetXCoordinates(), 
			     velocityGrid->GetYCoordinates(), 
			     velocityGrid->GetZCoordinates()};
  int nodeRes[3] = {coords[0]->GetNumberOfTuples(),
		    coords[1]->GetNumberOfTuples(),
		    coords[2]->GetNumberOfTuples()};
  int cellRes[3] = {nodeRes[0]-1,nodeRes[1]-1,nodeRes[2]-1};

  const int r = 1;//std::pow(2,refinement);
  // grid refinement comes here...
  const int subone = 0;//(refinement > 0 ? 1 : 0);
  int refNodeRes[3];
  refNodeRes[0] = nodeRes[0]*r - subone;
  refNodeRes[1] = nodeRes[1]*r - subone;
  refNodeRes[2] = nodeRes[2]*r - subone;

  int extent[6];
  velocityGrid->GetExtent(extent);
  int imin = extent[0] > globalExtent[0] ? numGhostLevels : 0;
  int imax = extent[1] < globalExtent[1] ? cellRes[0]-numGhostLevels : cellRes[0];
  int jmin = extent[2] > globalExtent[2] ? numGhostLevels : 0;
  int jmax = extent[3] < globalExtent[3] ? cellRes[1]-numGhostLevels : cellRes[1];
  int kmin = extent[4] > globalExtent[4] ? numGhostLevels : 0;
  int kmax = extent[5] < globalExtent[5] ? cellRes[2]-numGhostLevels : cellRes[2];

  for (int k = kmin; k < kmax; ++k) {
    double z = coords[2]->GetComponent(k,0);
    for (int j = jmin; j < jmax; ++j) {
      double y = coords[1]->GetComponent(j,0);
      for (int i = imin; i < imax; ++i) {
	double x = coords[0]->GetComponent(i,0);
	points->InsertNextPoint(x,y,z);	
      }
    }
  }
}

void generateSeedPointsPLIC(vtkRectilinearGrid *vofGrid,
			    int refinement,
			    vtkPoints *points,
			    int globalExtent[6],
			    int numGhostLevels)
{
  int index;
  vtkDataArray *vofArray =
    vofGrid->GetCellData()->GetArray("Data", index);
  if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }
  
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
    vofGrid->GetCellData()->GetArray("Data", index);
  if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }
  
  // vtkDataArray *data =
  //   vofGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);
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
}

#if 1

//static float interpolateSca(vtkDataArray *vofField,
//			      const int* res, const int idxCell[3],
//			      const double bcoords[3])

float4 vofCorrector(const float4 pos1, const float4 displacement,
		    vtkDataArray *vofField, vtkDataArray *coords[3],
		    const int res[3], const int ijk[3], const double pcoords[3])
{
  // -----------------------------------------------------------
  // compute gradient
  float vof[6];
  int ijk2[3] = {ijk[0], ijk[1], ijk[2]};

  //left
  if (ijk[0] > 0) {
    ijk2[0] = ijk[0]-1;
  }
  vof[0] = interpolateSca(vofField, res, ijk2, pcoords);
  ijk2[0] = ijk[0];
  //right
  if (ijk[0] < res[0]-1) {
    ijk2[0] = ijk[0]+1;
  }
  vof[1] = interpolateSca(vofField, res, ijk2, pcoords);
  ijk2[0] = ijk[0];
  //bottom
  if (ijk[1] > 0) {
    ijk2[1] = ijk[1]-1;
  }
  vof[2] = interpolateSca(vofField, res, ijk2, pcoords);
  ijk2[1] = ijk[1];
  //top
  if (ijk[1] < res[1]-1) {
    ijk2[1] = ijk[1]+1;
  }
  vof[3] = interpolateSca(vofField, res, ijk2, pcoords);
  ijk2[1] = ijk[1];
  //back
  if (ijk[2] > 0) {
    ijk2[2] = ijk[2]-1;
  }
  vof[4] = interpolateSca(vofField, res, ijk2, pcoords);
  ijk2[2] = ijk[2];
  //front
  if (ijk[2] < res[2]-1) {
    ijk2[2] = ijk[2]+1;
  }
  vof[5] = interpolateSca(vofField, res, ijk2, pcoords);
  ijk2[2] = ijk[2];

  float dx[6] = {(coords[0]->GetComponent(ijk[0]+1,0)-coords[0]->GetComponent(ijk[0]-1,0))/2.0f,
		 (coords[0]->GetComponent(ijk[0]+2,0)-coords[0]->GetComponent(ijk[0],  0))/2.0f,
		 (coords[1]->GetComponent(ijk[1]+1,0)-coords[1]->GetComponent(ijk[1]-1,0))/2.0f,
		 (coords[1]->GetComponent(ijk[1]+2,0)-coords[1]->GetComponent(ijk[1],  0))/2.0f,
		 (coords[2]->GetComponent(ijk[2]+1,0)-coords[2]->GetComponent(ijk[2]-1,0))/2.0f,
		 (coords[2]->GetComponent(ijk[2]+2,0)-coords[2]->GetComponent(ijk[2],  0))/2.0f};
  
  float3 grad = make_float3((vof[1]-vof[0])/(dx[0]+dx[1]),
  			    (vof[3]-vof[2])/(dx[2]+dx[3]),
  			    (vof[5]-vof[4])/(dx[4]+dx[5]));
  float3 ngrad = make_float3(0.0f);
  if (length(grad) > 0.0) {
   ngrad = normalize(grad);
  }
  float3 dd = make_float3((dx[0]+dx[1])/2.0f,
			  (dx[2]+dx[3])/2.0f,
			  (dx[4]+dx[5])/2.0f);
  float4 disp = make_float4(ngrad*dd,0.0f);
  return pos1 + disp;
  
  // float4 grad = make_float4((vof[1]-vof[0])/(dx[0]+dx[1]),
  // 			    (vof[3]-vof[2])/(dx[2]+dx[3]),
  // 			    (vof[5]-vof[4])/(dx[4]+dx[5]),0.0f);
  // float mag = length(make_float3(grad));
  // if (mag > 0.0f) {
  //   float4 disp = grad/mag;
  //   std::cout << "disp " << disp.x << " " << disp.y << " " << disp.z << std::endl;
  //   return pos1 + disp;
  // }  
  // return pos1;
} 

// iterative, solved with fixed point method - Newton's method can be viewed as such
// https://en.wikipedia.org/wiki/Fixed-point_iteration
// https://en.wikipedia.org/wiki/Trapezoidal_rule_%28differential_equations%29
void advectParticles(vtkRectilinearGrid *vofGrid,
		     vtkRectilinearGrid *velocityGrid,
		     std::vector<float4> &particles,
		     std::vector<float4> &velocities,		     
		     const float deltaT)
{
  int nodeRes[3];
  vofGrid->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  int index;
  vtkDataArray *velocityArray1 = velocityGrid->GetCellData()->GetArray("Data", index);
    if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }

  // vtkDataArray *velocityArray1 = velocityGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
  vtkDataArray *vofArray1 = vofGrid->GetCellData()->GetArray("Data", index);
  if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }
  // vtkDataArray *vofArray1 = vofGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

  vtkDataArray *coords[3] = {vofGrid->GetXCoordinates(),
			     vofGrid->GetYCoordinates(),
			     vofGrid->GetZCoordinates()};

  std::vector<float4>::iterator itp = particles.begin();
  std::vector<float4>::iterator itv = velocities.begin();
  for (; itp != particles.end() && itv != velocities.end(); ++itp, ++itv) {

    if (itp->w == 0.0f) {
      continue;
    }
    double x[3];
    int ijk[3];
    double pcoords[3];
    const int maxNumIter = 20;

    float4 pos0 = *itp;
    float4 velocity0 = *itv;
    // initial guess - forward Euler
    float4 pos1 = pos0 + deltaT*velocity0;
    float4 velocity1;
    
    {
      x[0] = pos1.x;
      x[1] = pos1.y;
      x[2] = pos1.z;      
      velocityGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
      int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      float f = vofArray1->GetComponent(idx, 0);
      if (f <= g_emf0) {
    	pos1 = vofCorrector(pos1, velocity0*deltaT, vofArray1, coords, cellRes, ijk, pcoords);
      }      
    }
      
    for (int i = 0; i < maxNumIter; ++i) {

      x[0] = pos1.x;
      x[1] = pos1.y;
      x[2] = pos1.z;
      velocityGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
      velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
      
      float4 velocity = (velocity0 + velocity1)/2.0f;
      pos1 = pos0 + deltaT*velocity;	
    }
    *itp = pos1;

    x[0] = itp->x;
    x[1] = itp->y;
    x[2] = itp->z;
    int particleInsideGrid = vofGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
    velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
    *itv = velocity1;
        
    if (particleInsideGrid) {
      int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      float f = vofArray1->GetComponent(idx, 0);

      itp->w = f;
      
      // if (f <= g_emf0) {
      // 	itp->w = 0.0f;
      // }
    }
  }
}
#else
// iterative, solved with fixed point method - Newton's method can be viewed as such
// https://en.wikipedia.org/wiki/Fixed-point_iteration
// https://en.wikipedia.org/wiki/Trapezoidal_rule_%28differential_equations%29
void advectParticles(vtkRectilinearGrid *vofGrid,
		     vtkRectilinearGrid *velocityGrid,
		     std::vector<float4> &particles,
		     std::vector<float4> &velocities,		     
		     const float deltaT)
{
  int nodeRes[3];
  vofGrid->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  int index;
  vtkDataArray *velocityArray1 = velocityGrid->GetCellData()->GetArray("Data", index);
    if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }

  // vtkDataArray *velocityArray1 = velocityGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);
  vtkDataArray *vofArray1 = vofGrid->GetCellData()->GetArray("Data", index);
  if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }
  // vtkDataArray *vofArray1 = vofGrid->GetCellData()->GetAttribute(vtkDataSetAttributes::SCALARS);

  std::vector<float4>::iterator itp = particles.begin();
  std::vector<float4>::iterator itv = velocities.begin();
  for (; itp != particles.end() && itv != velocities.end(); ++itp, ++itv) {

    if (itp->w == 0.0f) {
      continue;
    }
    double x[3];
    int ijk[3];
    double pcoords[3];
    const int maxNumIter = 20;

    float4 pos0 = *itp;
    float4 velocity0 = *itv;
    // initial guess - forward Euler
    float4 pos1 = pos0 + deltaT*velocity0;
    float4 velocity1;
      
    for (int i = 0; i < maxNumIter; ++i) {

      x[0] = pos1.x;
      x[1] = pos1.y;
      x[2] = pos1.z;
      velocityGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
      velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
      float4 velocity = (velocity0 + velocity1)/2.0f;
      pos1 = pos0 + deltaT*velocity;	
    }      
    *itp = pos1;

    x[0] = itp->x;
    x[1] = itp->y;
    x[2] = itp->z;
    int particleInsideGrid = vofGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
    velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
    *itv = velocity1;
        
    if (particleInsideGrid) {
      int idx = ijk[0] + ijk[1]*cellRes[0] + ijk[2]*cellRes[0]*cellRes[1];
      float f = vofArray1->GetComponent(idx, 0);

      if (f <= g_emf0) {
	itp->w = 0.0f;
      }
    }
  }
}
#endif

void advectParticles(vtkRectilinearGrid *velocityGrid[2],
		     std::vector<float4> &particles,
		     std::vector<float4> &velocities,
		     const float deltaT)
{
  int index;
  vtkDataArray *velocityArray0 = velocityGrid[0]->GetCellData()->GetArray("Data", index);
  if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }
  vtkDataArray *velocityArray1 = velocityGrid[1]->GetCellData()->GetArray("Data", index);
    if (index == -1) {
    std::cout << __LINE__ << ": Array not found!" << std::endl;
  }

  int nodeRes[3];
  velocityGrid->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  std::vector<float4>::iterator itp = particles.begin();
  std::vector<float4>::iterator itv = velocities.begin();
  for (; itp != particles.end() && itv != velocities.end(); ++itp, ++itv) {

    double x[3];
    int ijk[3];
    double pcoords[3];
    const int maxNumIter = 20;

    float4 pos0 = *itp;
    float4 velocity0 = *itv;
    // initial guess - forward Euler
    float4 pos1 = pos0 + deltaT*velocity0;
    float4 velocity1;
      
    for (int i = 0; i < maxNumIter; ++i) {

      x[0] = pos1.x;
      x[1] = pos1.y;
      x[2] = pos1.z;
      velocityGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
      velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
      float4 velocity = (velocity0 + velocity1)/2.0f;
      pos1 = pos0 + deltaT*velocity;	
    }      
    *itp = pos1;

    x[0] = itp->x;
    x[1] = itp->y;
    x[2] = itp->z;
    int particleInsideGrid = velocityGrid->ComputeStructuredCoordinates(x, ijk, pcoords);
    velocity1 = make_float4(interpolateVec(velocityArray1, cellRes, ijk, pcoords),0.0f);
    *itv = velocity1;
  }
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

void calcLabelPoints(vtkFloatArray *labels,
		     const int labelsRange[2],
		     std::vector<std::vector<int>> &labelPoints)
{
  const int numPoints = labels->GetNumberOfTuples();
  const int numLabels = labelsRange[1] - labelsRange[0] + 1;
  for (int i = 0; i < numPoints; ++i) {

    int pointLabel = labels->GetComponent(i, 0);
    if (pointLabel < labelsRange[0] || pointLabel > labelsRange[1]) {
      continue;
    }    
    int pointLabelIdx = pointLabel - labelsRange[0];
    labelPoints[pointLabelIdx].push_back(i);
  }
}

void calcLabelBounds(vtkPoints *points,
		     vtkFloatArray *labels,
		     const int labelsRange[2],
		     vtkRectilinearGrid *grid,
		     std::vector<std::array<int,6>> &labelBounds)
{
  const int numLabels = labelsRange[1] - labelsRange[0] + 1;
  std::vector<std::array<float,6>> labelBoundsTmp(numLabels);
  for (int i = 0; i < numLabels; ++i) {
    labelBoundsTmp[i][0] = labelBoundsTmp[i][2] = labelBoundsTmp[i][4] = 
      std::numeric_limits<float>::max();
    labelBoundsTmp[i][1] = labelBoundsTmp[i][3] = labelBoundsTmp[i][5] = 
      -1.0f*std::numeric_limits<float>::max();
  }

  // compute physical bounds of each label
  const int numPoints = points->GetNumberOfPoints();
  for (int i = 0; i < numPoints; ++i) {
    
    int pointLabel = labels->GetComponent(i, 0);
    if (pointLabel < labelsRange[0] || pointLabel > labelsRange[1]) {
      continue;
    }    

    int pointLabelIdx = pointLabel - labelsRange[0];

    double p[3];
    points->GetPoint(i, p);
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

// void generateBoundaries(vtkPoints *points,
// 			vtkFloatArray *labels,
// 			vtkRectilinearGrid *grid,			
// 			vtkPolyData *boundaries,
// 			const int refinement)
// {
//   if (points->GetNumberOfPoints() == 0) {
//     return;
//   }

//   const int subone = (refinement > 0 ? 1 : 0);
  
//   double range[2];
//   labels->GetRange(range, 0);
//   const int numLabels = std::ceil(range[1] - range[0] + 1.0f);

//   // stores indices of all points with a given label
//   std::vector<std::vector<int>> labelPoints(numLabels);
//   calcLabelPoints(labels, labelPoints);
  
//   // stores the index bounds of particles with a given label
//   std::vector<std::array<int,6>> labelBounds(numLabels);
//   calcLabelBounds(points, labels, grid, labelBounds);
  
//   const float isoValue = 0.501f;
//   int vertexID = 0;
//   std::vector<unsigned int> indices(0);
//   std::vector<float4> vertices(0);
//   vtkDataArray *coords[3] = {grid->GetXCoordinates(), 
// 			     grid->GetYCoordinates(), 
// 			     grid->GetZCoordinates()};

//   std::vector<int> labelOffsets(numLabels+1,0);
  
//   for (int i = 0; i < labelPoints.size(); ++i) {

//     if (labelPoints[i].size() == 0) {
//       labelOffsets[i+1] = vertices.size();
//       continue;
//     }

//     int ijk0[3] = {labelBounds[i][0],labelBounds[i][2],labelBounds[i][4]};
//     int ijk1[3] = {labelBounds[i][1],labelBounds[i][3],labelBounds[i][5]};

//     // this is a node-based grid so +1 for each dimension
//     int subNodeRes[3] = {ijk1[0]-ijk0[0]+1+1,
// 			 ijk1[1]-ijk0[1]+1+1,
// 			 ijk1[2]-ijk0[2]+1+1};

//     const int r = std::pow(2,refinement);
//     // grid refinement comes here...
//     subNodeRes[0] = subNodeRes[0]*r - subone;
//     subNodeRes[1] = subNodeRes[1]*r - subone;
//     subNodeRes[2] = subNodeRes[2]*r - subone;

//     vtkFloatArray *subcoords[3];

//     for (int n = 0; n < 3; n++) {
//       subcoords[n] = vtkFloatArray::New();
//       subcoords[n]->SetNumberOfComponents(1);
//       subcoords[n]->SetNumberOfTuples(subNodeRes[n]);

//       float xprev = coords[n]->GetComponent(ijk0[n], 0);    
//       int ires = (subNodeRes[n] + subone)/r;

//       for (int j = 0; j < ires-1; ++j) {
// 	float x = coords[n]->GetComponent(ijk0[n]+j+1, 0);
// 	float dx = (x - xprev)/r;
	
// 	for (int k = 0; k < r; ++k) {
// 	  subcoords[n]->SetValue(j*r+k, xprev+k*dx);
// 	}
// 	xprev = x;
//       }
//       subcoords[n]->SetValue(subNodeRes[n]-1, 
// 			     coords[n]->GetComponent(ijk1[n]+1,0)); 
//     }
    
//     vtkRectilinearGrid *subGrid = vtkRectilinearGrid::New();    
//     subGrid->SetDimensions(subNodeRes[0],subNodeRes[1],subNodeRes[2]);

//     subGrid->SetXCoordinates(subcoords[0]);
//     subGrid->SetYCoordinates(subcoords[1]);
//     subGrid->SetZCoordinates(subcoords[2]);

//     const int numElements = subNodeRes[0]*subNodeRes[1]*subNodeRes[2];
//     std::vector<float> field(numElements, 0.0f);

//     for (int j = 0; j < labelPoints[i].size(); ++j) {

//       double x[3];
//       points->GetPoint(labelPoints[i][j], x);
//       int ijk[3];
//       double pcoords[3];    
//       subGrid->ComputeStructuredCoordinates(x, ijk, pcoords);

//       int ids[8] =
// 	{ijk[0]   +  ijk[1]*subNodeRes[0]    +  ijk[2]*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]+1 +  ijk[1]*subNodeRes[0]    +  ijk[2]*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]+1 + (ijk[1]+1)*subNodeRes[0] +  ijk[2]*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]   + (ijk[1]+1)*subNodeRes[0] +  ijk[2]*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]   +  ijk[1]*subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]+1 +  ijk[1]*subNodeRes[0]    + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]+1 + (ijk[1]+1)*subNodeRes[0] + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1],
// 	 ijk[0]   + (ijk[1]+1)*subNodeRes[0] + (ijk[2]+1)*subNodeRes[0]*subNodeRes[1]};

//       field[ids[0]] += (1.0f-pcoords[0])*(1.0f-pcoords[1])*(1.0f-pcoords[2]);
//       field[ids[1]] += (     pcoords[0])*(1.0f-pcoords[1])*(1.0f-pcoords[2]);
//       field[ids[2]] += (     pcoords[0])*(     pcoords[1])*(1.0f-pcoords[2]);
//       field[ids[3]] += (1.0f-pcoords[0])*(     pcoords[1])*(1.0f-pcoords[2]);
//       field[ids[4]] += (1.0f-pcoords[0])*(1.0f-pcoords[1])*(     pcoords[2]);
//       field[ids[5]] += (     pcoords[0])*(1.0f-pcoords[1])*(     pcoords[2]);
//       field[ids[6]] += (     pcoords[0])*(     pcoords[1])*(     pcoords[2]);
//       field[ids[7]] += (1.0f-pcoords[0])*(     pcoords[1])*(     pcoords[2]);
//     }

//     extractSurface(field.data(), subNodeRes, subcoords, 0.501f, indices, vertices, vertexID);    

//     subGrid->Delete();

//     labelOffsets[i+1] = vertices.size();
//   }
  
//   vtkPoints *outputPoints = vtkPoints::New();
//   outputPoints->SetNumberOfPoints(vertices.size());
//   for (int i = 0; i < vertices.size(); ++i) {
    
//     double p[3] = {vertices[i].x,
//   		   vertices[i].y,
//   		   vertices[i].z};
//     outputPoints->SetPoint(i, p);        
//   }

//   vtkIdTypeArray *cells = vtkIdTypeArray::New();
//   cells->SetNumberOfComponents(1);
//   cells->SetNumberOfTuples(indices.size()/3*4);
//   for (int i = 0; i < indices.size()/3; ++i) {
//     cells->SetValue(i*4+0,3);
//     cells->SetValue(i*4+1,indices[i*3+0]);
//     cells->SetValue(i*4+2,indices[i*3+1]);
//     cells->SetValue(i*4+3,indices[i*3+2]);
//   }

//   vtkCellArray *outputTriangles = vtkCellArray::New();
//   outputTriangles->SetNumberOfCells(indices.size()/3);
//   outputTriangles->SetCells(indices.size()/3, cells);

//   vtkShortArray *boundaryLabels = vtkShortArray::New();
//   boundaryLabels->SetName("Labels");
//   boundaryLabels->SetNumberOfComponents(1);
//   boundaryLabels->SetNumberOfTuples(outputPoints->GetNumberOfPoints());

//   for (int i = 0; i < numLabels; ++i) {
//     for (int j = labelOffsets[i]; j < labelOffsets[i+1]; ++j) {
//       boundaryLabels->SetValue(j, i);
//     }
//   }

//   boundaries->SetPoints(outputPoints);
//   boundaries->SetPolys(outputTriangles);
//   boundaries->GetPointData()->AddArray(boundaryLabels);

// }

void resamplePointsOnGrid(const std::vector<int> &labelPoints,
			  vtkPoints *points,
			  vtkRectilinearGrid *subGrid,
			  const int subNodeRes[3],
			  std::vector<float> &field)
{
  for (int j = 0; j < labelPoints.size(); ++j) {

    double x[3];
    points->GetPoint(labelPoints[j], x);
    int ijk[3];
    double pcoords[3];    
    subGrid->ComputeStructuredCoordinates(x, ijk, pcoords);

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

void generateCoords(vtkDataArray *coords, const int subNodeRes,
		    const int subone, const int r,
		    const int ijk0, const int ijk1,
		    vtkFloatArray **subcoords)
{
  (*subcoords) = vtkFloatArray::New();
  (*subcoords)->SetNumberOfComponents(1);
  (*subcoords)->SetNumberOfTuples(subNodeRes);

  float xprev = coords->GetComponent(ijk0, 0);    
  int ires = (subNodeRes + subone)/r;

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

void generateBoundaries(vtkPoints *points,
			vtkPoints *neighborPoints,
			vtkFloatArray *labels,
			vtkFloatArray *neighborLabels,
			vtkRectilinearGrid *grid,			
			vtkPolyData *boundaries,
			const int refinement)
{
  if (points->GetNumberOfPoints() == 0) {
    return;
  }
  const int subone = (refinement > 0 ? 1 : 0);
  
  double range[2];
  labels->GetRange(range, 0);
  const int numLabels = std::ceil(range[1] - range[0] + 1.0f); // buggy if range[0] > 0!!!
  const int labelsRange[2] = {range[0], range[1]};

  // stores indices of all points with a given label
  std::vector<std::vector<int>> labelPoints(numLabels);
  std::vector<std::vector<int>> labelNeighborPoints(numLabels);

  for (int i = 0; i < numLabels; ++i) {
    labelPoints[i].clear();
    labelNeighborPoints[i].clear();
  }
  calcLabelPoints(labels, labelsRange, labelPoints);
  if (neighborPoints != nullptr && neighborPoints->GetNumberOfPoints() > 0) {
    calcLabelPoints(neighborLabels, labelsRange, labelNeighborPoints);
  }
  
  // stores the index bounds of particles with a given label
  std::vector<std::array<int,6>> labelBounds(numLabels);
  std::vector<std::array<int,6>> labelNeighborBounds(numLabels);
  
  calcLabelBounds(points, labels, labelsRange, grid, labelBounds);
  if (neighborPoints != nullptr && neighborPoints->GetNumberOfPoints() > 0) {
    calcLabelBounds(neighborPoints, neighborLabels, labelsRange, grid, labelNeighborBounds);
  }
  else {
    labelNeighborBounds = labelBounds;
  }

  const float isoValue = 0.501f;
  int vertexID = 0;
  std::vector<unsigned int> indices(0);
  std::vector<float4> vertices(0);
  vtkDataArray *coords[3] = {grid->GetXCoordinates(), 
			     grid->GetYCoordinates(), 
			     grid->GetZCoordinates()};

  std::vector<int> labelOffsets(numLabels+1,0);
  
  for (int i = 0; i < labelPoints.size(); ++i) {

    if (labelPoints[i].size() == 0) {
      labelOffsets[i+1] = vertices.size();
      continue;
    }

    int ijk0[3] = {std::min(labelBounds[i][0], labelNeighborBounds[i][0]),
		   std::min(labelBounds[i][2], labelNeighborBounds[i][2]),
		   std::min(labelBounds[i][4], labelNeighborBounds[i][4])};
    int ijk1[3] = {std::max(labelBounds[i][1], labelNeighborBounds[i][1]),
		   std::max(labelBounds[i][3], labelNeighborBounds[i][3]),
		   std::max(labelBounds[i][5], labelNeighborBounds[i][5])};

    // this is a node-based grid so +1 for each dimension
    int subNodeRes[3] = {ijk1[0]-ijk0[0]+1+1,
			 ijk1[1]-ijk0[1]+1+1,
			 ijk1[2]-ijk0[2]+1+1};    

    const int r = std::pow(2,refinement);
    // grid refinement comes here...
    subNodeRes[0] = subNodeRes[0]*r - subone;
    subNodeRes[1] = subNodeRes[1]*r - subone;
    subNodeRes[2] = subNodeRes[2]*r - subone;
    
    int extent[6] = {0, subNodeRes[0]-1,
		     0, subNodeRes[1]-1,
		     0, subNodeRes[2]-1};

    // --
    if (labelBounds[i][0] > labelNeighborBounds[i][0]) extent[0] += r;
    if (labelBounds[i][1] < labelNeighborBounds[i][1]) extent[1] -= r;
    if (labelBounds[i][2] > labelNeighborBounds[i][2]) extent[2] += r;
    if (labelBounds[i][3] < labelNeighborBounds[i][3]) extent[3] -= r;
    if (labelBounds[i][4] > labelNeighborBounds[i][4]) extent[4] += r;
    if (labelBounds[i][5] < labelNeighborBounds[i][5]) extent[5] -= r;    

    vtkFloatArray *subcoords[3];
    for (int n = 0; n < 3; n++) {
      generateCoords(coords[n], subNodeRes[n], subone, r, ijk0[n],
      		     ijk1[n], &subcoords[n]);
    }
    
    vtkRectilinearGrid *subGrid = vtkRectilinearGrid::New();    
    subGrid->SetDimensions(subNodeRes[0],subNodeRes[1],subNodeRes[2]);

    subGrid->SetXCoordinates(subcoords[0]);
    subGrid->SetYCoordinates(subcoords[1]);
    subGrid->SetZCoordinates(subcoords[2]);

    const int numElements = subNodeRes[0]*subNodeRes[1]*subNodeRes[2];
    std::vector<float> field(numElements, 0.0f);

    resamplePointsOnGrid(labelPoints[i], points, subGrid, subNodeRes, field);    
    resamplePointsOnGrid(labelNeighborPoints[i], neighborPoints, subGrid, subNodeRes, field);

    extractSurface(field.data(), subNodeRes, subcoords, extent, 0.501f, indices, vertices, vertexID);    

    subGrid->Delete();

    labelOffsets[i+1] = vertices.size();
  }
  
  vtkPoints *outputPoints = vtkPoints::New();
  outputPoints->SetNumberOfPoints(vertices.size());
  for (int i = 0; i < vertices.size(); ++i) {
    
    double p[3] = {vertices[i].x,
  		   vertices[i].y,
  		   vertices[i].z};
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

  vtkShortArray *boundaryLabels = vtkShortArray::New();
  boundaryLabels->SetName("Labels");
  boundaryLabels->SetNumberOfComponents(1);
  boundaryLabels->SetNumberOfTuples(outputPoints->GetNumberOfPoints());

  for (int i = 0; i < numLabels; ++i) {
    for (int j = labelOffsets[i]; j < labelOffsets[i+1]; ++j) {
      boundaryLabels->SetValue(j, i);
    }
  }

  boundaries->SetPoints(outputPoints);
  boundaries->SetPolys(outputTriangles);
  boundaries->GetPointData()->AddArray(boundaryLabels);

}
