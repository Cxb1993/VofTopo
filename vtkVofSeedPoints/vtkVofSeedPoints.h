#ifndef __vtkVofSeedPoints_h
#define __vtkVofSeedPoints_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkExecutive.h"

class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;
class vtkPoints;
class vtkShortArray;

#include "vtkMath.h"

class vtkVofSeedPoints : public vtkPolyDataAlgorithm
{
 public:
  vtkSetMacro(Refinement, int);
  vtkGetMacro(Refinement, int);

  vtkSetMacro(Reseed, int);
  vtkGetMacro(Reseed, int);

  static vtkVofSeedPoints *New();
  vtkTypeMacro(vtkVofSeedPoints, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  void AddSourceConnection(vtkAlgorithmOutput* input);
  void RemoveAllSources();

 protected:
  vtkVofSeedPoints();
  ~vtkVofSeedPoints();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); 

 private:

  bool DataOnCells;
  int Refinement;
  vtkPoints *OutputSeeds;
  vtkIntArray *Connectivity;
  vtkShortArray *Coords;
  int Reseed;
};
#endif
