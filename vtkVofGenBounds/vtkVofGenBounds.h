#ifndef __vtkVofGenBounds_h
#define __vtkVofGenBounds_h

#include "vtkPolyDataAlgorithm.h" //superclass

class vtkInformation;
class vtkInformationVector;

#include "vtkMath.h"

class vtkVofGenBounds : public vtkPolyDataAlgorithm
{
 public:

  static vtkVofGenBounds *New();
  vtkTypeMacro(vtkVofGenBounds, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  void AddSourceConnection(vtkAlgorithmOutput* input);
  void RemoveAllSources();

 protected:
  vtkVofGenBounds();
  ~vtkVofGenBounds();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); 

 private:
};
#endif
