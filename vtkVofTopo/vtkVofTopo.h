#ifndef __vtkVofTopo_h
#define __vtkVofTopo_h

#include "vtkMultiBlockDataSetAlgorithm.h"
#include <map>
#include <vector>
#include <string>
#include "vofTopology.h"

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

  vtkGetMacro(Refinement, int);
  vtkSetMacro(Refinement, int);

  vtkGetMacro(ParticleStoreFreq, int);
  vtkSetMacro(ParticleStoreFreq, int);

  vtkGetMacro(StoreIntermParticles, int);
  vtkSetMacro(StoreIntermParticles, int);

  vtkGetMacro(StoreIntermBoundaries, int);
  vtkSetMacro(StoreIntermBoundaries, int);

  vtkGetMacro(IntegrationMethod, int);
  vtkSetMacro(IntegrationMethod, int);

  vtkGetMacro(PLICCorrection, int);
  vtkSetMacro(PLICCorrection, int);

  vtkGetMacro(VOFCorrection, int);
  vtkSetMacro(VOFCorrection, int);

  vtkGetMacro(RK4NumSteps, int);
  vtkSetMacro(RK4NumSteps, int);  

  vtkGetMacro(EMF0, double);
  vtkSetMacro(EMF0, double);  

  vtkGetMacro(EMF1, double);
  vtkSetMacro(EMF1, double);  

  vtkGetMacro(SeedByPLIC, int);
  vtkSetMacro(SeedByPLIC, int);  

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
  void InitParticles(vtkRectilinearGrid *vof, vtkPolyData *seeds);
  void AdvectParticles(vtkRectilinearGrid *vof[2],
		       vtkRectilinearGrid *velocity[2]);
  void ExchangeParticles();
  void ExtractComponents(vtkRectilinearGrid *vof,
			 vtkRectilinearGrid *components);

  void LabelAdvectedParticles(vtkRectilinearGrid *components,
			      std::vector<float> &labels);
  void TransferParticleDataToSeeds(std::vector<float> &particleData,
				   const std::string arrayName,
				   vtkPolyData *dst);
  void TransferIntermParticlesToSeeds(std::vector<std::vector<float4>> &particles,
				      std::vector<std::vector<int>> &ids,
				      std::vector<std::vector<short>> &procs);
  
  void GenerateBoundaries(vtkPolyData *boundaries, 
			  vtkPolyData *boundarySeeds,
			  vtkPolyData *seeds);

  void ExchangeBoundarySeedPoints(vtkPolyData *boundarySeeds, vtkPolyData *seeds);

  void GenerateIntParticles(vtkPolyData *intParticles);
  void GenerateIntBoundaries(vtkPolyData *intermBoundaries);

  void CreateScalarField(vtkRectilinearGrid *grid,
			 const std::vector<float4> &particles,
			 vtkRectilinearGrid *scalarField);

  std::vector<double> InputTimeValues;
  
  int InitTimeStep; // time t0
  int TargetTimeStep; // time t1 = t0+T

  int TimestepT0;
  int TimestepT1;
  
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
  int LocalExtent[NUM_SIDES];
  int LocalExtentNoGhosts[NUM_SIDES];

  // Seeds
  int Refinement;
  vtkPolyData *Seeds;
  bool SeedPointsProvided;

  // Particles
  std::vector<float4> Particles;
  std::vector<int> ParticleIds;
  std::vector<short> ParticleProcs;
  std::vector<float> Uncertainty;

  // Intermediate particles
  std::vector<std::vector<float4>> IntermParticles;
  std::vector<std::vector<int>> IntermParticleIds;
  std::vector<std::vector<short>> IntermParticleProcs;
  std::vector<float> IntermParticlesTimeStamps;
  
  // Temporal boundaries
  vtkPolyData *Boundaries;
  // Intermediate boundaries
  std::vector<std::vector<int>> IntermBoundaryIndices;
  std::vector<std::vector<float4>> IntermBoundaryVertices;
  std::vector<std::vector<float4>> IntermBoundaryNormals;
  std::vector<std::vector<int>> PrevLabelPoints;
  
  int VertexID;
  
  // Caching  
  bool UseCache;
  int LastLoadedTimestep;

  // Vof and velocity
  vtkRectilinearGrid *VofGrid[2];
  vtkRectilinearGrid *VelocityGrid[2];

  // Misc
  int IntegrationMethod; // 0 - Heun, 1 - RK4
  int PLICCorrection;
  int VOFCorrection;
  int RK4NumSteps;
  int StoreIntermParticles;
  int StoreIntermBoundaries;

  double EMF0;
  double EMF1;

  int SeedByPLIC;

  int ParticleStoreFreq;
};

#endif
