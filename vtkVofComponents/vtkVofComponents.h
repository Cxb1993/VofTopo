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
  virtual int RequestInformation(vtkInformation* request,
				 vtkInformationVector** inputVector,
				 vtkInformationVector* outputVector);

  virtual int RequestData(vtkInformation *, 
			  vtkInformationVector **, 
			  vtkInformationVector *);

 private:
  vtkMPIController *Controller;
};
#endif
