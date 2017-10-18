/*********************
* Date of creation 17/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
*
**************************************/


#ifndef _MOVING_ON_H
#define _MOVING_ON_H

#include <stdlib.h>

void moving_on_routine(int numberOfParticles, double timeStep, double temperature, double *diffusionMatrix, double *additionalForces, double *stochasticDisplacement, double *generalisedCoordinates);

#endif // _MOVING_ON_H
