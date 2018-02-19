#ifndef _FORCES_H
#define _FORCES_H

#include <stdlib.h>

#include "particles.h"

void force_torque_summation(double *additionalForces,double *generalisedCoordinates, int numberOfCells, int *forceList, int numberOfForces,	environmentVariables conditions);


#endif //_FORCES_H
