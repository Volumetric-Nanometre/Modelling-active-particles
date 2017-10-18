
/*********************
* Date of creation 17/10/2017
* Author:
* Contact:
* Other Authors: N/A
**************************************
* History
**************************************/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "stochastic_force.h"


void stochastic_displacement_creation(int numberOfParticles, double *stochasticDisplacement)
{
    memset(stochasticDisplacement,'0', 6 * numberOfParticles );
}
