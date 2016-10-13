#ifndef CORRESPONDENCE_H
#define CORRESPONDENCE_H

#include "vtkRectilinearGrid.h"
#include "vtkPolyData.h"

void GetCorrespondences(vtkRectilinearGrid *sourceGrid,
			vtkRectilinearGrid *destinationGrid,
			vtkRectilinearGrid *correspondences,
			const float distThreshold = 0.0f);

#endif//CORRESPONDENCE_H
