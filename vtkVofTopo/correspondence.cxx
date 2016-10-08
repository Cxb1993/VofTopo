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

//----------------------------------------------------------------------------
void GetCorrespondences(vtkRectilinearGrid *sourceGrid,
			vtkRectilinearGrid *destinationGrid,
			vtkRectilinearGrid *correspondences)

{
  std::cout << "bangla" << std::endl;
  typedef uint2 assignment_e_t;
  typedef std::tuple<intervol::AssignmentsPlain<assignment_e_t>,
		     std::vector<float3>,
		     std::vector<double>
		     > rslt_t;
  const size_t maxValue = 255;

  int res[3];
  sourceGrid->GetDimensions(res);
  res[0] -= 1;
  res[1] -= 1;
  res[2] -= 1;
  
  int3 volDim0 = make_int3(res[0],res[1],res[2]);
  int3 volDim1 = volDim0;
  
  intervol::DataPosSupply<float3, uint32_t> dataPosSupply0;
  intervol::DataPosSupply<float3, uint32_t> dataPosSupply1;

  int numCells = res[0]*res[1]*res[2];
  
  std::vector<float> data0(numCells,0.0f);
  std::vector<float> data1(numCells,0.0f);

  vtkFloatArray *sourceArray = vtkFloatArray::SafeDownCast(sourceGrid->GetCellData()->GetArray(0));
  vtkFloatArray *destinationArray = vtkFloatArray::SafeDownCast(destinationGrid->GetCellData()->GetArray(0));

  for (int i = 0; i < numCells; ++i) {
    data0[i] = sourceArray->GetComponent(i,0);
    data1[i] = destinationArray->GetComponent(i,0);
  }
  
  // loadData(data0, data1, volDim0, volDim1, po, maxValue);

  // const float3 ext0 = getDomainSize(po.configFNames.front());
  // const float3 ext1 = getDomainSize(po.configFNames.back());
  // const float3 voxelSize0 = make_float3(1.0f);//getVoxelSize(po.configFNames.front());
  // const float3 voxelSize1 = make_float3(1.0f);//getVoxelSize(po.configFNames.back());
  
  // const double normFac = 1.0f;//getNormFac(ext0, ext1);   
  
  const float3 f0 = make_float3(1.0f);// voxelSize0*normFac;
  const float3 f1 = make_float3(1.0f);// voxelSize1*normFac;
  
  intervol::initPosSupply(dataPosSupply0.pos, dataPosSupply0.supply,
  			  data0, volDim0, f0);
  intervol::initPosSupply(dataPosSupply1.pos, dataPosSupply1.supply,
  			  data1, volDim1, f1);

  
  double termDist = 0.0;
  intervol::TermDistTime term(termDist, 15);
  intervol::IterLog log;
  intervol::Plan_random plan;
  
  rslt_t rslt = intervol::run<intervol::AssignmentsPlain<assignment_e_t>>(dataPosSupply0,
									  dataPosSupply1,
									  term, log, plan);
  
  auto assignments = std::get<0>(rslt);
  std::vector<float3> pos = std::get<1>(rslt);
  std::vector<double> sth = std::get<2>(rslt);
  
  for(size_t i = 0; i< assignments.size(); i++) {

    std::vector<std::pair<float,int>> sortedWeights;
    sortedWeights.clear();

    for(auto e_it = assignments.cbegin(i); e_it != assignments.cend(i); e_it++) {
      sortedWeights.push_back(std::pair<float,int>(e_it->y,e_it->x));
    }

    std::sort(sortedWeights.begin(), sortedWeights.end(), sortByFirst);
    std::pair<float,int> corr = sortedWeights.back();

    if (i < 100) {
      float3 p = pos[corr.second];
      std::cout << "p = " << p.x << " " << p.y << " " << p.z << std::endl;
    }
    
    // int idx = sortedWeights.size()-1;

    // float dist2 = length(pos[i]-pos[sortedWeights[idx].second]);
    // std::cout << dist2 << std::endl;
    // while (idx >= 0 && dist2 > 0.012) {
    //   --idx;
    //   if (idx < 0) break;
    //   dist2 = length(pos[i]-pos[sortedWeights[idx].second]);
    // }

    // if (idx >= 0) {
    //   idlist[1] = sortedWeights[idx].second;
    //   lines->InsertNextCell(2,idlist);
    //   weights->InsertNextValue(sortedWeights[idx].first);
    // }

    // idlist[1] = sortedWeights.back().second;
    // lines->InsertNextCell(2,idlist);
    // weights->InsertNextValue(sortedWeights.back().first);
  }

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
