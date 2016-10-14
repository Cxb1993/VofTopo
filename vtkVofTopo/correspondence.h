#ifndef CORRESPONDENCE_H
#define CORRESPONDENCE_H

#include "vtkRectilinearGrid.h"
#include "vtkPolyData.h"

void GetCorrespondences(vtkRectilinearGrid *sourceGrid,
			vtkRectilinearGrid *destinationGrid,
			vtkRectilinearGrid *correspondences,
			vtkPolyData *corrs);

#endif//CORRESPONDENCE_H
