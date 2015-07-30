#include "vofTopology.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include <iostream>
#include <map>
#include <vector>
#include <limits>

namespace
{
  typedef struct {
    int x;
    int y;
    int z;
  } int3_t;

  inline int3_t operator/(int3_t a, int b)
  {
    int3_t c = {a.x/b, a.y/b, a.z/b};
    return c;
  }

  bool int3_t_compare(const int3_t &a, const int3_t &b)
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
		  std::map<int3_t, int, bool(*)(const int3_t &a, const int3_t &b)> &seedPos,
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
	    int3_t pos = {cell_x*subdiv + xr, 
			  cell_y*subdiv + yr, 
			  cell_z*subdiv + zr};
	    seedPos[pos] = seedIdx;
	    ++seedIdx;

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

  std::map<int3_t, int, bool(*)(const int3_t &a, const int3_t &b)> seedPos(int3_t_compare);
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
    connectivity->SetTuple3(seedIdx, conn[0], conn[1], conn[2]);
    coords->SetTuple3(seedIdx, x, y, z);
  }
}


void advectParticles(vtkRectilinearGrid *inputVof,
		     vtkRectilinearGrid *inputVelocity,
		     std::vector<float4> &particles,
		     const float deltaT)
{
  int nodeRes[3];
  inputVof->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  vtkDataArray *velocityArray = inputVelocity->GetCellData()->GetAttribute(vtkDataSetAttributes::VECTORS);

  std::vector<float4>::iterator it;
  for (it = particles.begin(); it != particles.end(); ++it) {
    
    double x[3] = {it->x, it->y, it->z};
    int ijk[3];
    double pcoords[3];
    int particleInsideGrid = inputVof->ComputeStructuredCoordinates(x, ijk, pcoords);

    if (particleInsideGrid) {
      float4 velocity = make_float4(interpolateVec(velocityArray, cellRes, ijk, pcoords), 0.0f);
      *it = *it + velocity*deltaT;
    }
  }
}

// multiprocess
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
