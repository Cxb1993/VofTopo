#ifndef __vtkVofAdvect_h
#define __vtkVofAdvect_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared
#include "vtkPolyDataAlgorithm.h"
#include "vtkExecutive.h"
#include "vtkMath.h"

#include <vector>
#include <set>
#include <map>

#include "helper_math.h"
#include "voftopo.h"

class vtkPolyData;
class vtkCellArray;
class vtkFloatArray;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;
class vtkRectilinearGrid;
class vtkMPIController;

class vtkVofAdvect : public vtkPolyDataAlgorithm
{
 public:
  static vtkVofAdvect *New();
  vtkTypeMacro(vtkVofAdvect, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  void AddSourceConnection(vtkAlgorithmOutput* input);
  void RemoveAllSources();

  // GUI -------------------------------
  vtkGetMacro(StartTimeStep, int);
  vtkSetMacro(StartTimeStep, int);
  //~GUI -------------------------------

 protected:
  vtkVofAdvect();
  ~vtkVofAdvect();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  int RequestInformation(vtkInformation* request,
			 vtkInformationVector** inputVector,
			 vtkInformationVector* outputVector);

  int RequestUpdateExtent(vtkInformation* request,
			  vtkInformationVector** inputVector,
			  vtkInformationVector* outputVector);
  
  // Generate output
  int RequestData(vtkInformation *, 
		  vtkInformationVector **, 
		  vtkInformationVector *);

  int ProcessData(vtkRectilinearGrid *inputVofGrid,
		  vtkRectilinearGrid *inputVelocityGrid);

 private:

  std::vector<double> InputTimeValues;
  double StartTime;
  double TerminationTime;
  int StartTimeStep;
  int TerminationTimeStep;
  int CurrentTimeStep;
  bool FirstIteration;

  std::vector<f3u1_t> Particles;
  std::vector<short> PProcessId;

  vtkDataArray *XCoords;
  vtkDataArray *YCoords;
  vtkDataArray *ZCoords;

  // multiprocess
  vtkMPIController* Controller;
  double MyBounds[6];
  double GlobalBounds[6];
  std::vector<std::vector<int> > NeighborProcesses;
  int NumNeighbors;
  int GlobalExtents[6];
};
#endif
