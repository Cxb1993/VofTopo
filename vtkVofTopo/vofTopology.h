#ifndef VOFTOPOLOGY_H
#define VOFTOPOLOGY_H

#include "vtkRectilinearGrid.h"
#include "vtkPoints.h"
#include "vtkShortArray.h"
#include "vtkIntArray.h"
#include "vtkMPIController.h"
#include <vector>
#include "helper_math.h"

void generateSeedPoints(vtkRectilinearGrid *input,
			int refinement,
			vtkPoints *points,
			vtkIntArray *connectivity,
			vtkShortArray *coords);

void advectParticles(vtkRectilinearGrid *inputVof,
		     vtkRectilinearGrid *inputVelocity,
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

  
#endif//VOFTOPOLOGY_H
