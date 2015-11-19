#ifndef MARCHINGCUBES_CPU_H
#define MARCHINGCUBES_CPU_H

#include "vtkDataArray.h"
#include <vector_types.h>
#include <vector>

void extractSurface(const float* volume, 
		    const unsigned* resolution,
		    vtkDataArray *coords[3],
		    const int gridOffset[3],
		    const float isoValue,		    	    
		    std::vector<unsigned int>& indices,
		    std::vector<float4>& vertices,
		    int &vertexID);

#endif//MARCHINGCUBES_CPU_H
