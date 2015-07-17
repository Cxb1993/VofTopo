#ifndef __vtkVofTopo_h
#define __vtkVofTopo_h

#include "vtkPolyDataAlgorithm.h"
#include <vector>

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
  //~GUI -------------------------------


protected:
  vtkVofTopo();
  ~vtkVofTopo();

  int FillInputPortInformation(int, vtkInformation*);
  int RequestInformation(vtkInformation *,
			 vtkInformationVector **,
			 vtkInformationVector *);
  int RequestUpdateExtent(vtkInformation *,
			  vtkInformationVector **,
			  vtkInformationVector *);
  int RequestData(vtkInformation *,
		  vtkInformationVector **,
		  vtkInformationVector *);

private:

  vtkVofTopo(const vtkVofTopo&);  // Not implemented.
  void operator=(const vtkVofTopo&);  // Not implemented.

  std::vector<double> InputTimeValues;
  
  int InitTimeStep; // time t0
  int TargetTimeStep; // time t1 = t0+T
  int CurrentTimeStep;
  // we can iterate over t0 or t1
  enum IterationType {IterateInit, IterateTarget};
  IterationType IterType;
  bool FirstIteration;

  // for data sets without or with incorrect time stamp information
  double TimeStepDelta;
};

#endif
