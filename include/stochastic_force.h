/*********************
* Date of creation 17/10/2017
* Author: Oliver Hinds
* Contact:
* Other Authors: N/A
**************************************
* History
*
**************************************/

#ifndef _STOCHASTIC_FORCE_H
#define _STOCHASTIC_FORCE_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>


void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, long int tSeed);

#endif // _STOCHASTIC_FORCE_H
