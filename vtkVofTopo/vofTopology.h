#ifndef VOFTOPOLOGY_H
#define VOFTOPOLOGY_H

#include "vtkRectilinearGrid.h"
#include "vtkPoints.h"
#include "vtkShortArray.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkMPIController.h"
#include <vector>
#include <algorithm>
#include <map>
#include "helper_math.h"

int findClosestTimeStep(double requestedTimeValue,
			const std::vector<double>& timeSteps);


void generateSeedPoints(vtkRectilinearGrid *input,
			int refinement,
			vtkPoints *points,
			vtkIntArray *connectivity,
			vtkShortArray *coords,
			int onlyOnInterface);

void generateSeedPointsPLIC(vtkRectilinearGrid *input,
			    int refinement,
			    vtkPoints *points,
			    vtkIntArray *connectivity,
			    vtkShortArray *coords,
			    int onlyOnInterface);

void advectParticles(vtkRectilinearGrid *inputVof[2],
		     vtkRectilinearGrid *inputVelocity[2],
		     std::vector<float4> &particles,
		     const float deltaT);

// multiprocess
void findGlobalExtents(std::vector<int> &allExtents, 
		       int globalExtents[6]);
void findGlobalBounds(std::vector<double> &allBounds, 
		      double globalBounds[6]);
void findNeighbors(const int myExtents[6], 
		   const int globalExtents[6], 
		   const std::vector<int> &allExtents,
		   std::vector<std::vector<int> > &neighbors);

int outOfBounds(const float4 particle, const double bounds[6], const double globalBounds[6]);
int withinBounds(const float4 particle, const double bounds[6]);
  
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

// components
static const double g_emf0 = 0.000001;
static const double g_emf1 = 0.999999;
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

void prepareLabelsToSend(std::vector<std::vector<int> > &NeighborProcesses,
			 const int myExtent[6], int cellRes[3], vtkFloatArray *labels,
			 std::vector<std::vector<float4> > &labelsToSend);

void unifyLabelsInProcess(std::vector<std::vector<int> > &NeighborProcesses,
			  const int myExtent[6], int cellRes[3], vtkFloatArray *labels,
			  std::vector<std::vector<float4> > &labelsToRecv,
			  std::vector<int> &labelOffsets, int processId,
			  std::vector<int> &allLabels);

void unifyLabelsInDomain(std::vector<int> &allLabelUnions, int numAllLabels,
			 std::vector<int> &allLabels, vtkFloatArray *labels,
			 std::vector<int> &labelOffsets, int processId);

void generateBoundaries(vtkPoints *points,
			vtkFloatArray *labels,
			vtkIntArray *connectivity,
			vtkShortArray *coords,
			vtkFloatArray *velos,
			vtkPolyData *boundaries);

void regenerateBoundaries(vtkPoints *points,
			  vtkFloatArray *labels,
			  vtkIntArray *connectivity,
			  vtkShortArray *coords,
			  int currentTimeStep,			  
			  std::vector<float3> &vertices,
			  std::vector<float3> &ivertices,
			  std::vector<int> &indices,
			  std::vector<int> &splitTimes);

void regenerateBoundaries(vtkPoints *points, vtkFloatArray *labels,
			  vtkIntArray *connectivity, vtkShortArray *coords,
			  int currentTimeStep, vtkPolyData *boundaries);

void smoothSurface(std::vector<float3>& vertices,
		   std::vector<int>& indices);
void smoothSurface(vtkPoints *vertices,
		   vtkCellArray *cells);

void mergePatches(std::vector<float3>& vertices,
		  std::vector<float3>& ivertices,
		  std::vector<int>& indices);

void mergePatches(vtkPolyData *boundaries);

void computeParticleVelocities(std::vector<float4> &particles,
			       vtkRectilinearGrid *velocityGrid,
			       std::vector<float4> &velocities);

#endif//VOFTOPOLOGY_H
