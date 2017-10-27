/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: Oliver Hinds
**************************************
* Change History
**************************************/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <CL\cl.h>

#include "particles.h"
#include "diffusionmatrix.h"
#include "stochastic_force.h"
#include "moving_on.h"


double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;

int gDebug = 0;

int main(int argc, char *argv[])
{
    //
    // Check Debug mode
    //



    if(argc >= 2)
    {
        gDebug = 1;
        printf("Warning debug mode entered. Press any key to continue...\n" );
        getchar();
    }

    FILE *output = fopen("output.txt","w");

    if(output == NULL)
    {
        printf("-Error %d : %s\n", errno, strerror( errno ) );
        return -errno;
    }



    int numberOfParticles = 0;



    particleVariables* particles = NULL;
    //
    // Call function to read in particle data
    //
    if( ( numberOfParticles = particle_read_in( &particles ) ) <= 0)
    {
        getchar();
        return numberOfParticles;
    }

    printf("Data read in success\n" );

    //
    // Create generalised coordinates
    //

    double *generalisedCoordinates = NULL ;

    if( (generalisedCoordinates = generalised_coordinate_creation( numberOfParticles, particles) ) == NULL )
    {
        free( particles );
        particles = NULL ;
        getchar();
        return errno;
    }

    //---------------------------- DEBUG------------------------------//
    //
    // Prints the generalisedCoordinates to a file for inspection
    //
    if( gDebug == 1 && generalisedCoordinates != NULL)
    {
        FILE *genCoordOutput = fopen( "genCoord_output.txt","w");

        for(int i = 0; i < 6 * numberOfParticles; i++)
        {
            fprintf(genCoordOutput, "%e\n", generalisedCoordinates[i]);
        }
        fclose (genCoordOutput);
    }
    //---------------------------END---------------------------------//


    //
    // Allocate memory required for the program.
    // Requires: Diffusion matrix, stochastic displacement,
    //          additional forces
    //




    double *diffusionMatrix = NULL ;
    double *stochasticWeighting = NULL;
    double *additionalForces = NULL ;
    double *stochasticDisplacement = NULL;

    if( (diffusionMatrix = calloc( pow( 6 * numberOfParticles, 2), sizeof *diffusionMatrix) ) == NULL)
    {
        free( particles );
        particles = NULL ;
        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;
        printf("-Error %d : %s\n", errno, strerror( errno ) );

        getchar();
        return -errno;
    }

    if( (stochasticWeighting = calloc( pow( 6 * numberOfParticles, 2), sizeof *stochasticWeighting) ) == NULL)
    {
        free( particles );
        particles = NULL ;
        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;
        free( diffusionMatrix );
        diffusionMatrix = NULL;
        printf("-Error %d : %s\n", errno, strerror( errno ) );

        getchar();
        return -errno;
    }

    if( (stochasticDisplacement = calloc( 6 * numberOfParticles, sizeof *stochasticDisplacement) ) == NULL)
    {
        free( particles );
        particles = NULL ;
        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;
        free( diffusionMatrix );
        diffusionMatrix = NULL;
        free( stochasticWeighting );
        stochasticWeighting = NULL;
        printf("-Error %d : %s\n", errno, strerror( errno ) );

        getchar();
        return -errno;
    }

    if( (additionalForces = calloc( 6 * numberOfParticles, sizeof *additionalForces) ) == NULL)
    {
        free( particles );
        particles = NULL ;
        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;
        free( diffusionMatrix );
        diffusionMatrix = NULL;
        free( stochasticDisplacement );
        stochasticDisplacement = NULL;
        free( stochasticWeighting );
        stochasticWeighting = NULL;
        printf("-Error %d : %s\n", errno, strerror( errno ) );

        getchar();
        return -errno;
    }


    double temperature = 298; // K
    double viscosity = 8.9E-4; //N m^-2 s
    double radius = 1E-6; // m
    double currentTime = 0;
    double deltaTime = 1E-3; // Seconds
    double endTime = 1E-2; // Seconds

    while(currentTime<=endTime)
    {
        //
        // Create diffusion matrix
        //

        diffusion_matrix_creation( numberOfParticles, diffusionMatrix, stochasticWeighting, particles, temperature, viscosity, radius);

        //---------------------------- DEBUG------------------------------//
        //
        // Prints the diffusionMatrix to a file for inspection
        //
        if( gDebug == 1 && diffusionMatrix != NULL)
        {
            currentTime = endTime+1;
            FILE *matrixOutput = fopen( "matrix_output.txt","w");

            for(int i = 0; i < 6 * numberOfParticles; i++)
            {
                for(int j = 0; j < 6 * numberOfParticles; j++)
                {
                    fprintf(matrixOutput, "%e\t", diffusionMatrix[i * 6 * numberOfParticles + j]);
                }
                fprintf(matrixOutput, "\n");

            }
            fclose (matrixOutput);
        }
        //---------------------------END---------------------------------//


        //
        // Create the stochastic displacement
        //

		time_t tSeed;
		time(&tSeed);

        stochastic_displacement_creation( numberOfParticles, stochasticWeighting, stochasticDisplacement, tSeed );

        //
        // Include additional forces
        //


        //
        // Calculate time step.
        //

        moving_on_routine(numberOfParticles, deltaTime, temperature, diffusionMatrix, additionalForces, stochasticDisplacement, generalisedCoordinates);
        fprintf(output, "%lf\t", currentTime);
        for(int i = 0; i < 6 * numberOfParticles; i++)
        {
            fprintf(output, "%e\t", generalisedCoordinates[i]);
        }
        fprintf(output, "\n");

        currentTime+=deltaTime; // time step
    }

    //
    // Free memory
    //

    fclose (output);

    if( particles != NULL )
    {
        free( particles );
        particles = NULL;
    }
    if( diffusionMatrix != NULL )
    {
        free( diffusionMatrix );
        particles = NULL;
    }
    if( generalisedCoordinates != NULL )
    {
        free( generalisedCoordinates );
        particles = NULL;
    }
    if( stochasticWeighting != NULL )
    {
        free( stochasticWeighting );
        particles = NULL;
    }
    if( stochasticDisplacement != NULL )
    {
        free( stochasticDisplacement );
        particles = NULL;
    }

    return 0;
}
