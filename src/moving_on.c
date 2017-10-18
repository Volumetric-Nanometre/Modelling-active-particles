
/*********************
* Date of creation 17/10/2017
* Author: MichaelO'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
**************************************/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "moving_on.h"

extern double gBoltzmannConst;

void moving_on_routine(int numberOfParticles, double timeStep, double temperature, double *diffusionMatrix, double *additionalForces, double *stochasticDisplacement, double *generalisedCoordinates)
{
    double frontConst = timeStep / ( temperature * gBoltzmannConst);

    for(int i = 0; i < 6 * numberOfParticles; i++)
    {
        //
        // Calculate  D F * deltat
        //
        double temp = 0;

        //#pragma omp parallel for reduction(+ : temp)
        for(int j = 0; j < 6 * numberOfParticles; j++)
        {
            temp += frontConst*diffusionMatrix[i * numberOfParticles + j] * additionalForces[j];
        }

        //
        // Add the diffusion displacement to the current position and the
        // random stochastic displacement to get the final displacement.
        //
        // Replace the original generalisedCoordinates with the new coordinates
        //

        generalisedCoordinates[i] = temp + stochasticDisplacement[i] + generalisedCoordinates[i];
    }

    //
    // Calculate dD/dy *deltat ( currently this term is ignored )
    //

}
