#ifndef CORRESPONDENCE_H
#define CORRESPONDENCE_H

#include "vtkRectilinearGrid.h"

void GetCorrespondences(vtkRectilinearGrid *sourceGrid,
			vtkRectilinearGrid *destinationGrid,
			vtkRectilinearGrid *correspondences);

#endif//CORRESPONDENCE_H
