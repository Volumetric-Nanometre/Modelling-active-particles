/*********************
* Date of creation 17/10/2017
* Author: Oliver Hinds
* Contact:
**************************************
* History
*
**************************************/

#ifndef _STOCHASTIC_FORCE_H
#define _STOCHASTIC_FORCE_H

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>


void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, gsl_rng *tSeed, double timestep);

#endif // _STOCHASTIC_FORCE_H
