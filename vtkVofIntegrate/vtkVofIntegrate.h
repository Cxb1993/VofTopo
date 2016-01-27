#ifndef __vtkVofIntegrate_h
#define __vtkVofIntegrate_h

//#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "helper_math.h"
#include <map>
#include <vector>

class vtkMPIController;
class vtkRectilinearGrid;
class vtkPolyData;
class vtkFloatArray;

class VTK_EXPORT vtkVofIntegrate : public vtkPolyDataAlgorithm //vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkVofIntegrate* New();
  vtkTypeMacro(vtkVofIntegrate, vtkPolyDataAlgorithm); //vtkMultiBlockDataSetAlgorithm);
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
  //~GUI -------------------------------

protected:
  vtkVofIntegrate();
  ~vtkVofIntegrate();

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

  vtkVofIntegrate(const vtkVofIntegrate&);  // Not implemented.
  void operator=(const vtkVofIntegrate&);  // Not implemented.

  void GetGlobalContext(vtkInformation *inInfo);
  void InitParticles(vtkRectilinearGrid *vof);
  void InitVelocities(vtkRectilinearGrid *velocity);
  void AdvectParticles(vtkRectilinearGrid *vof[2],
		       vtkRectilinearGrid *velocity[2]);
  void ExchangeParticles();
  void ExtractComponents(vtkRectilinearGrid *vof,
			 vtkRectilinearGrid *components);

  void LabelAdvectedParticles(vtkRectilinearGrid *components,
			      std::vector<float> &labels);
  void TransferLabelsToSeeds(std::vector<float> &particleLabels);

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

  double Incr;

  // for data sets without or with incorrect time stamp information
  double TimeStepDelta;

  // Multiprocess
  vtkMPIController* Controller;
  static const int NUM_SIDES = 6;
  double LocalBounds[NUM_SIDES];
  double GlobalBounds[NUM_SIDES];
  double BoundsNoGhosts[NUM_SIDES];
  std::vector<std::vector<int> > NeighborProcesses;
  int NumNeighbors;
  int NumGhostLevels;
  int GlobalExtent[NUM_SIDES];

  // Seeds
  int Refinement;
  vtkPolyData *Seeds;

  // Particles
  std::vector<float4> Particles;
  std::vector<float4> Velocities;
  std::vector<int> ParticleIds;
  std::vector<short> ParticleProcs;
  
  // Caching  
  bool UseCache;
  int LastLoadedTimestep;

  // Vof and velocity
  vtkRectilinearGrid *VofGrid[2];
  vtkRectilinearGrid *VelocityGrid[2];
};

#endif
