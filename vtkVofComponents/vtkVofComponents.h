#ifndef __vtkVofComponents_h
#define __vtkVofComponents_h

#include "vtkRectilinearGridAlgorithm.h" //superclass

#include <vector>

class vtkMPIController;

class vtkVofComponents : public vtkRectilinearGridAlgorithm
{
 public:
  static vtkVofComponents *New();
  vtkTypeMacro(vtkVofComponents, vtkRectilinearGridAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  void AddSourceConnection(vtkAlgorithmOutput* input);
  void RemoveAllSources();

 protected:
  vtkVofComponents();
  ~vtkVofComponents();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  // Generate output
  int RequestUpdateExtent(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);
  virtual int RequestData(vtkInformation *, 
			  vtkInformationVector **, 
			  vtkInformationVector *);

 private:
  vtkMPIController *Controller;
  static const int NUM_SIDES = 6;
  double LocalBounds[NUM_SIDES];
  double GlobalBounds[NUM_SIDES];
  std::vector<std::vector<int> > NeighborProcesses;
  int NumNeighbors;
  int NumGhostLevels;
  int GlobalExtent[NUM_SIDES];
  
  void GetGlobalContext(vtkInformation *inInfo);
  
  void ExtractComponents(vtkRectilinearGrid *vof,
			 vtkRectilinearGrid *components);
};

#endif
