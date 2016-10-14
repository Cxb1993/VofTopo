#include "correspondence.h"

#include "vtkFloatArray.h"
#include "vtkCellData.h"

#include <algorithm>
#include <utility>

//----------------------------------------------------------------------------
#define USE_OMP
#define NO_FLANN

#include <boost/predef.h>
#include "volData/vec.h"

#include "volData/VolDataHandlerCUDA.h"
#include "volData/equalizeSum.h"
#ifdef USE_OMP
#include <omp.h>
#endif//USE_OMP
#include "intervol_init.h"
#include "intervol_eval.h"
#include "intervol_io.h"
#include "intervol_app_common.h"
#include "intervol.h"
#include "intervol_po.h"
#include "intervol_plan.h"
#include "intervol_Log.h"
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
static bool sortByFirst(std::pair<float,int> a, std::pair<float,int> b)
{
  return a.first < b.first;
}

// //----------------------------------------------------------------------------
// void GetCorrespondences(vtkRectilinearGrid *sourceGrid,
// 			vtkRectilinearGrid *destinationGrid,
// 			vtkRectilinearGrid *correspondences,
// 			vtkPolyData *corrs,
// 			const float distThreshold)

// {
//   typedef uint2 assignment_e_t;
//   typedef std::tuple<intervol::AssignmentsPlain<assignment_e_t>,
// 		     std::vector<float3>,
// 		     std::vector<double>
// 		     > rslt_t;

//   int cellRes[3];
//   sourceGrid->GetDimensions(cellRes);
//   cellRes[0] -= 1;
//   cellRes[1] -= 1;
//   cellRes[2] -= 1;
  
//   int3 volDim0 = make_int3(cellRes[0],cellRes[1],cellRes[2]);
//   int3 volDim1 = volDim0;
  
//   int numCells = cellRes[0]*cellRes[1]*cellRes[2];

//   const int maxValue = 1 << 16;
  
//   std::vector<int> data0(numCells,0.0f);
//   std::vector<int> data1(numCells,0.0f);

//   vtkFloatArray *sourceArray = vtkFloatArray::SafeDownCast(sourceGrid->GetCellData()->GetArray(0));
//   vtkFloatArray *destinationArray = vtkFloatArray::SafeDownCast(destinationGrid->GetCellData()->GetArray(0));

//   for (int i = 0; i < numCells; ++i) {
//     data0[i] = sourceArray->GetComponent(i,0)*maxValue;
//   }
//   for (int i = 0; i < numCells; ++i) {
//     data1[i] = destinationArray->GetComponent(i,0)*maxValue;
//   }

  
//   const float cellScale = 1.0f/cellRes[0];
//   const float rcellScale = cellRes[0];
  
//   const float3 f0 = make_float3(cellScale);
//   const float3 f1 = make_float3(cellScale);
  
//   intervol::DataPosSupply<float3, uint32_t> dataPosSupply0;
//   intervol::DataPosSupply<float3, uint32_t> dataPosSupply1;
  
//   intervol::initPosSupply(dataPosSupply0.pos, dataPosSupply0.supply,
//   			  data0, volDim0, f0);
//   intervol::initPosSupply(dataPosSupply1.pos, dataPosSupply1.supply,
//   			  data1, volDim1, f1);

//   intervol::TermBasic term(1024);
//   intervol::IterLog log;
//   intervol::Plan_random plan;
  
//   rslt_t rslt = intervol::run<intervol::AssignmentsPlain<assignment_e_t>>(dataPosSupply0,
// 									  dataPosSupply1,
// 									  term, log, plan);
  
//   auto assignments = std::get<0>(rslt);
//   std::vector<float3> pos = std::get<1>(rslt);
//   std::vector<double> sth = std::get<2>(rslt);

//   vtkFloatArray *correspondenceVectors = vtkFloatArray::New();
//   correspondenceVectors->SetName("CorrespondenceVectors");
//   correspondenceVectors->SetNumberOfComponents(3);
//   correspondenceVectors->SetNumberOfTuples(numCells);
//   correspondenceVectors->FillComponent(0,0.0f);
//   correspondenceVectors->FillComponent(1,0.0f);
//   correspondenceVectors->FillComponent(2,0.0f);

  
//   for(size_t i = 0; i< assignments.size(); i++) {
    
//     // idlist[0] = i;
//     float3 p0 = pos[i];
//     int3 src = make_int3(pos[i].x*rcellScale,
// 			 pos[i].y*rcellScale,
// 			 pos[i].z*rcellScale);

//     float3 v = make_float3(0.0f);
//     float totalWeight = 0.0f;
    
//     for(auto e_it = assignments.cbegin(i); e_it != assignments.cend(i); e_it++) {

//       float3 p1 = pos[e_it->x];
//       float3 vec = p1 - p0;

//       if (length(vec) > distThreshold) {
// 	continue;
//       }
      
//       float dist2 = dot(vec,vec);
      
//       float weight = 1.0f;
//       if (dist2 > 0.0f) {
//       	weight = 1.0f/dist2;
//       }
//       else {
//       	weight = 1.0f/(0.125f*cellScale*cellScale);
//       }
//       weight = weight*std::sqrt(e_it->y);

//       v += vec*weight;
//       totalWeight += weight;
      
//       // break;
//     }
//     if (totalWeight > 0.0f) {
//       v = v/totalWeight;
//     }

//     int idx_src = int(src.x) + int(src.y)*cellRes[0] + int(src.z)*cellRes[0]*cellRes[1];
//     correspondenceVectors->SetTuple3(idx_src, v.x, v.y, v.z);
//   }

//   correspondences->GetCellData()->AddArray(correspondenceVectors);

//   assignments.free();
// }

void GetCorrespondences(vtkRectilinearGrid *sourceGrid,
			vtkRectilinearGrid *destinationGrid,
			vtkRectilinearGrid *correspondences,
			vtkPolyData *corrs)

{
  typedef uint2 assignment_e_t;
  typedef std::tuple<intervol::AssignmentsPlain<assignment_e_t>,
		     std::vector<float3>,
		     std::vector<double>
		     > rslt_t;

  int cellRes[3];
  sourceGrid->GetDimensions(cellRes);
  cellRes[0] -= 1;
  cellRes[1] -= 1;
  cellRes[2] -= 1;
  
  int3 volDim0 = make_int3(cellRes[0],cellRes[1],cellRes[2]);
  int3 volDim1 = volDim0;
  
  int numCells = cellRes[0]*cellRes[1]*cellRes[2];

  const int maxValue = 1 << 4;
  
  std::vector<int> data0(numCells,0.0f);
  std::vector<int> data1(numCells,0.0f);

  vtkFloatArray *sourceArray = vtkFloatArray::SafeDownCast(sourceGrid->GetCellData()->GetArray(0));
  vtkFloatArray *destinationArray = vtkFloatArray::SafeDownCast(destinationGrid->GetCellData()->GetArray(0));

  for (int i = 0; i < numCells; ++i) {
    data0[i] = sourceArray->GetComponent(i,0)*maxValue;
  }
  for (int i = 0; i < numCells; ++i) {
    data1[i] = destinationArray->GetComponent(i,0)*maxValue;
  }  
  
  const float3 f0 = make_float3(1.0f/cellRes[0]);
  const float3 f1 = make_float3(1.0f/cellRes[0]);

  intervol::DataPosSupply<float3, uint32_t> dataPosSupply0;
  intervol::DataPosSupply<float3, uint32_t> dataPosSupply1;

  intervol::initPosSupply(dataPosSupply0.pos, dataPosSupply0.supply,
  			  data0, volDim0, f0);
  intervol::initPosSupply(dataPosSupply1.pos, dataPosSupply1.supply,
  			  data1, volDim1, f1);

  double termDist=0.;
  intervol::TermDistTime term(termDist, 10e8);
  intervol::IterLog<intervol::IterLogElem> log;
  intervol::Plan_random plan;

  rslt_t rslt = intervol::run<intervol::AssignmentsPlain<assignment_e_t>>(dataPosSupply0,
									  dataPosSupply1,
									  term, log, plan);  
  auto assignments = std::get<0>(rslt);
  std::vector<float3> pos = std::get<1>(rslt);
  std::vector<double> sth = std::get<2>(rslt);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i=0; i< pos.size(); i++) {
    double p[3] = {pos[i].x, pos[i].y, pos[i].z};
    points->InsertNextPoint(p);
  }

  corrs->SetPoints(points);

  vtkIdType idlist[2];
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
  weights->SetNumberOfComponents(1);
  weights->SetName("weights");

  
  for(size_t i = 0; i< assignments.size(); i++) {

    idlist[0] = i;
    
    for(auto e_it = assignments.cbegin(i); e_it != assignments.cend(i); e_it++) {

      idlist[1] = e_it->x;
      lines->InsertNextCell(2,idlist);
      weights->InsertNextValue(e_it->y);
    }
  }

  corrs->SetLines(lines);
  corrs->GetCellData()->AddArray(weights);

  assignments.free();
}


// reductionFactors
// selectSliceParams
// configFNames
// 0       /home/karchgz/dev/intervol/configs/oil_funs0017_resampled_raw.config
// 1       /home/karchgz/dev/intervol/configs/oil_funs0017_raw.config
// transFuncFNames
// timeSteps
// 0       0
// termTime = 1e+08
// termImprovement = 1e-06
// termIterations = 18446744073709551615
// distanceLimit = 1e-06
// preNormData = 0
// outDir = out_interp
// alignCenters = 0
// nIntermedates = 10
// termDistRef = 
// doWriteResult = 1
