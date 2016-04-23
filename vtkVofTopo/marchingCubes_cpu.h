#ifndef MARCHINGCUBES_CPU_H
#define MARCHINGCUBES_CPU_H

#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include <vector>

#include "vectors_cuda.h"

void extractSurface(const float* volume, 
		    const int* resolution,
		    vtkFloatArray *coords[3],
		    const int extent[6],
		    const float isoValue,		    	    
		    std::vector<int>& indices,
		    std::vector<float4>& vertices,
		    int &vertexID);

#endif//MARCHINGCUBES_CPU_H
