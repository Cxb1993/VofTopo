#ifndef __vtkVofTopo_h
#define __vtkVofTopo_h

#include "vtkPolyDataAlgorithm.h"
#include "helper_math.h"
#include <vector>

class vtkMPIController;
class vtkRectilinearGrid;
class vtkPolyData;
class vtkFloatArray;

class VTK_EXPORT vtkVofTopo : public vtkPolyDataAlgorithm
{
public:
  static vtkVofTopo* New();
  vtkTypeMacro(vtkVofTopo, vtkPolyDataAlgorithm);
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
  void GenerateNewBoundaries(vtkPolyData *boundaries);


  std::vector<double> InputTimeValues;
  
  int InitTimeStep; // time t0
  int TargetTimeStep; // time t1 = t0+T
  int CurrentTimeStep;
  
  // we can iterate over t0 or t1
  static const int IterateOverInit = 0;
  static const int IterateOverTarget = 1;
  int IterType;
  bool FirstIteration;

  // Visualization type
  static const int LabelComponents = 0;
  static const int LabelSplitTime = 1;
  int LabelType;

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

  // Caching  
  bool UseCache;
  int LastComputedTimeStep;
};

#endif
