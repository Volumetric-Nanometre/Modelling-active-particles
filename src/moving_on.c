
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
    //printf("front const %e\n", frontConst);

    for(int i = 0; i < 6 * numberOfParticles; i++)
    {

        double temp = 0;
        //
        // Calculate  D F * deltat
        //
        //#pragma omp parallel for reduction(+ : temp)
        for(int j = 0; j < 6 * numberOfParticles; j++)
        {
            temp += frontConst*diffusionMatrix[i * 6*numberOfParticles + j] * additionalForces[j];
    //        printf("temp %e\t diff %e \tforce : %e     i %d j %d\n",temp,diffusionMatrix[i * 6*numberOfParticles + j],additionalForces[j], i,j );
        //    getchar();
        }
        //
        // Add the diffusion displacement to the current position and the
        // random stochastic displacement to get the final displacement.
        //
        temp += stochasticDisplacement[i];
        //
        // Calculate the velocity of the particle
        //

        //printf("final temp %d = %e\n",i,temp );
    //    getchar();
        //velocities[i] = temp /conditions->deltaTime;
        //
        // Replace the original generalisedCoordinates with the new coordinates
        //
        generalisedCoordinates[i] = temp + generalisedCoordinates[i];
    }

    //
    // Calculate dD/dy *deltat ( currently this term is ignored )
    //

}
