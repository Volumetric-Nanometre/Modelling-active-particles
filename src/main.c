/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* Change History
**************************************/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "particles.h"
#include "diffusionmatrix.h"
#include "stochastic_force.h"


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


    int numberOfParticles = 0;

    int deltaTime = 1E-3; // Seconds

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




    double *diffusionMatrix = NULL ;

    double temperature = 273; // Kelvin
    double viscosity = 1;
    double radius = 1;

    //
    // Create diffusion matrix
    //

    if( (diffusionMatrix = diffusion_matrix_creation( numberOfParticles, particles ,temperature, viscosity, radius)) == NULL )
    {

        free( particles );
        particles = NULL ;
        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;
        getchar();
        return errno;
    }

    printf("Diffusion matrix created\n" );

    //---------------------------- DEBUG------------------------------//
    //
    // Prints the diffusionMatrix to a file for inspection
    //
    if( gDebug == 1 && diffusionMatrix != NULL)
    {
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
    // Create the stochastic force
    //

    double *stochasticForce = NULL;

    if( (stochasticForce = stochastic_force_creation( numberOfParticles ) ) == NULL )
    {
        free( particles );
        particles = NULL ;

        free( diffusionMatrix );
        diffusionMatrix = NULL ;

        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;

        getchar();
        return errno;
    }



    //
    // Include additional forces
    //

    double *additionalForces = NULL;

    if( (additionalForces = calloc(6*numberOfParticles, sizeof *additionalForces) ) == NULL )
    {
        free( particles );
        particles = NULL ;

        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;

        free( diffusionMatrix );
        diffusionMatrix = NULL ;

        free( stochasticForce );
        stochasticForce = NULL ;

        getchar();
        return errno;
    }

    //
    // Calculate time step.
    //






    getchar();

    //
    // Free memory
    //

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
    if( stochasticForce != NULL )
    {
        free( stochasticForce );
        particles = NULL;
    }

    return 0;
}
