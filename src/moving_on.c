
/*********************
* Date of creation 17/10/2017
* Author: MichaelO'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
**************************************/

#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <omp.h>

#include "moving_on.h"
#include "particles.h"

extern double gBoltzmannConst;

void moving_on_routine(int numberOfParticles, environmentVariables *conditions, double *diffusionMatrix, double *additionalForces, double *stochasticDisplacement, double *generalisedCoordinates, double *velocities)
{
    double frontConst = conditions->deltaTime / ( conditions->temperature * gBoltzmannConst);

    for(int i = 0; i < 6 * numberOfParticles; i++)
    {

        double temp = 0;
        //
        // Calculate  D F * deltat
        //
        #pragma omp parallel for reduction(+ : temp)
        for(int j = 0; j < 6 * numberOfParticles; j++)
        {
            temp += frontConst*diffusionMatrix[i * numberOfParticles + j] * additionalForces[j];
        }
        //
        // Add the diffusion displacement to the current position and the
        // random stochastic displacement to get the final displacement.
        //
        temp += stochasticDisplacement[i] + generalisedCoordinates[i];
        //
        // Calculate the velocity of the particle
        //
        velocities[i] = (temp - generalisedCoordinates[i]) * conditions->deltaTime;
        //
        // Replace the original generalisedCoordinates with the new coordinates
        //
        generalisedCoordinates[i] = temp;
    }

    //
    // Calculate dD/dy *deltat ( currently this term is ignored )
    //

}
