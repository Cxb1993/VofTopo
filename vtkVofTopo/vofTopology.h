#ifndef VOFTOPOLOGY_H
#define VOFTOPOLOGY_H

#include "vtkRectilinearGrid.h"
#include "vtkPoints.h"
#include "vtkShortArray.h"
#include "vtkIntArray.h"
#include <vector>
#include "helper_math.h"

void generateSeedPoints(vtkRectilinearGrid *input,
			int refinement,
			vtkPoints *points,
			vtkIntArray *connectivity,
			vtkShortArray *coords);

void advectParticles(vtkRectilinearGrid *inputVof,
		     vtkRectilinearGrid *inputVelocity,
		     std::vector<float4> &particles,
		     const float deltaT);

#endif//VOFTOPOLOGY_H
