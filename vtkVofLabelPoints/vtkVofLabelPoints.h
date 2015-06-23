#ifndef __vtkVofLabelPoints_h
#define __vtkVofLabelPoints_h

#include "vtkMultiBlockDataSetAlgorithm.h" //superclass

class vtkMPIController;

class vtkVofLabelPoints : public vtkMultiBlockDataSetAlgorithm
{
 public:
  static vtkVofLabelPoints *New();
  vtkTypeMacro(vtkVofLabelPoints, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  void AddSourceConnection(vtkAlgorithmOutput* input);
  void RemoveAllSources();

 protected:
  vtkVofLabelPoints();
  ~vtkVofLabelPoints();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  // Generate output
  int RequestData(vtkInformation *, 
		  vtkInformationVector **, 
		  vtkInformationVector *);

 private:
  vtkMPIController* Controller;
};
#endif
