#ifndef _INITIAL_VELOCITIES_H
#define _INITIAL_VELOCITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <omp.h>

#include "particles.h"

void initialVelocities(int numberOfParticles, particleVariables *particles, environmentVariables *conditions, time_t tSeed);


#endif
