#ifndef __vtkVofTopo_h
#define __vtkVofTopo_h

// #include "vtkPolyDataAlgorithm.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "helper_math.h"
#include <map>
#include <vector>

class vtkMPIController;
class vtkRectilinearGrid;
class vtkPolyData;
class vtkFloatArray;

class VTK_EXPORT vtkVofTopo : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkVofTopo* New();
  vtkTypeMacro(vtkVofTopo, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // GUI -------------------------------
  vtkGetMacro(InitTimeStep, int);
  vtkSetMacro(InitTimeStep, int);

  vtkGetMacro(TargetTimeStep, int);
  vtkSetMacro(TargetTimeStep, int);

  vtkGetMacro(TimeStepDelta, double);
  vtkSetMacro(TimeStepDelta, double);

  vtkGetMacro(IterType, int);
  vtkSetMacro(IterType, int);

  vtkGetMacro(Refinement, int);
  vtkSetMacro(Refinement, int);

  vtkGetMacro(ComputeComponentLabels, int);
  vtkSetMacro(ComputeComponentLabels, int);

  vtkGetMacro(ComputeSplitTime, int);
  vtkSetMacro(ComputeSplitTime, int);

  vtkGetMacro(StepIncr, int);
  vtkSetMacro(StepIncr, int);

//~GUI -------------------------------

protected:
  vtkVofTopo();
  ~vtkVofTopo();

  int FillInputPortInformation(int, vtkInformation*);
  int RequestInformation(vtkInformation*,
			 vtkInformationVector**,
			 vtkInformationVector*);
  int RequestUpdateExtent(vtkInformation*,
			  vtkInformationVector**,
			  vtkInformationVector*);
  int RequestData(vtkInformation*,
		  vtkInformationVector**,
		  vtkInformationVector*);

private:

  vtkVofTopo(const vtkVofTopo&);  // Not implemented.
  void operator=(const vtkVofTopo&);  // Not implemented.

  void GetGlobalContext(vtkInformation *inInfo);
  void GenerateSeeds(vtkRectilinearGrid *vof);
  void InitParticles();
  void AdvectParticles(vtkRectilinearGrid *vof,
		       vtkRectilinearGrid *velocity);
  void ExchangeParticles();
  void ExtractComponents(vtkRectilinearGrid *vof,
			 vtkRectilinearGrid *components);
  void LabelAdvectedParticles(vtkRectilinearGrid *components,
			      std::vector<float> &labels);
  void TransferLabelsToSeeds(std::vector<float> &particleLabels);

  void GenerateBoundaries(vtkPolyData *boundaries);
  void GenerateTemporalBoundaries(vtkPolyData *boundaries, bool lastStep);


  std::vector<double> InputTimeValues;
  
  int InitTimeStep; // time t0
  int TargetTimeStep; // time t1 = t0+T

  int TimestepT0;
  int TimestepT1;
  
  // we can iterate over t0 or t1
  static const int ITERATE_OVER_INIT = 0;
  static const int ITERATE_OVER_TARGET = 1;
  int IterType;

  // Visualization type
  int ComputeComponentLabels;
  int ComputeSplitTime;

  // Integration direction
  static const int INTEGRATION_FORWARD = 1;
  static const int INTEGRATION_BACKWARD = -1;
  int StepIncr;

  // for data sets without or with incorrect time stamp information
  double TimeStepDelta;

  // Multiprocess
  vtkMPIController* Controller;
  static const int NUM_SIDES = 6;
  double LocalBounds[NUM_SIDES];
  double GlobalBounds[NUM_SIDES];
  std::vector<std::vector<int> > NeighborProcesses;
  int NumNeighbors;

  // Seeds
  int Refinement;
  vtkPolyData *Seeds;

  // Particles
  std::vector<float4> Particles;
  std::vector<int> ParticleIds;
  std::vector<short> ParticleProcs;
  
  // Temporal boundaries
  vtkPolyData *Boundaries;

  // Caching  
  bool UseCache;
  int LastLoadedTimestep;

  // Vof and velocity
  vtkRectilinearGrid *VofGrid[2];
  vtkRectilinearGrid *VelocityGrid[2];
};

#endif
