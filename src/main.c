/*********************
* Date of creation 09/10/2017
* Authors: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* Change History
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "particles.h"
#include "diffusionmatrix.h"
#include "stochastic_force.h"
#include "moving_on.h"
<<<<<<< HEAD
=======
//#include "initial_velocities.h"
>>>>>>> aa2f782f24ba8d13debbe77abecf6a3c4c54ab03
#include "forces.h"

double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;
double gGrav = 9.80665; // m s^-2

int gDebug = 0;

int main(int argc, char *argv[])
{
    //
    // Check Debug mode
    //



    if(argc > 1)
    {
        gDebug = 1;
        printf("Warning debug mode entered. Press any key to continue...\n" );
        getchar();
    }

    FILE *output = fopen("../bin/output.txt","w");

    if(output == NULL)
    {
        printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
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
        FILE *genCoordOutput = fopen("../bin/genCoord_output.txt","w");

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
    //          additional forces, stochasticWeighting,
    //          velocities
    //


    double *diffusionMatrix = NULL ;
    double *stochasticWeighting = NULL;
    double *additionalForces = NULL;
    double *stochasticDisplacement = NULL;

    diffusionMatrix = calloc( pow( 6 * numberOfParticles, 2), sizeof *diffusionMatrix) ;
    stochasticWeighting = calloc( pow( 6 * numberOfParticles, 2), sizeof *stochasticWeighting);
    stochasticDisplacement = calloc( 6 * numberOfParticles, sizeof *stochasticDisplacement);
    additionalForces = calloc( 6 * numberOfParticles, sizeof *additionalForces);

    if(  diffusionMatrix==NULL  || stochasticWeighting==NULL || stochasticDisplacement==NULL || additionalForces==NULL)
    {
        free( particles );
        particles = NULL ;

        free( generalisedCoordinates );
        generalisedCoordinates = NULL ;

        if( diffusionMatrix!=NULL)
        free( diffusionMatrix );
        diffusionMatrix = NULL;

        if( stochasticDisplacement!=NULL)
        free( stochasticDisplacement );
        stochasticDisplacement = NULL;

        if( stochasticWeighting!=NULL)
        free( stochasticWeighting );
        stochasticWeighting = NULL;

        if( additionalForces!=NULL)
        free( additionalForces );
        additionalForces = NULL;


        printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);

        return -errno;
    }


    //
    // Allocate the environmental conditions and nano particle
    // characteristics
    //

	environmentVariables conditions;
    conditions.temperature = 298; // K
    conditions.viscosity = 8.9E-4; //N m^-2 s
    conditions.radius = 1E-6; // m
    conditions.currentTime = 0;
    conditions.deltaTime = 1E-5; // Seconds
    conditions.endTime = 2; // Seconds
	conditions.mass = (4/3) * gPi * pow(conditions.radius,3)*19320; // kg - density of gold

    //
    //  Choose forces to be included
    //

    int numberOfForces = 2; // must be at least 1, with the force none chosen

    int forceList[2] = {0,0}; // repulsive and van der waals



	time_t tSeed1;
	time(&tSeed1);
	long int tSeed = -1*(long int) tSeed1;

    //
    // Loop through time, output each time step to a file.
    //
    while(conditions.currentTime<=conditions.endTime)
    {
        //
        // Create diffusion matrix
        //

        diffusion_matrix_creation( numberOfParticles, diffusionMatrix, stochasticWeighting, generalisedCoordinates, &conditions);

        //---------------------------- DEBUG------------------------------//
        //
        // Prints the diffusionMatrix to a file for inspection
        //
        if( gDebug == 1 && diffusionMatrix != NULL)
        {
            //conditions.currentTime = conditions.endTime+1;
            FILE *matrixOutput = fopen("../bin/matrix_output.txt","a");

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

        stochastic_displacement_creation( numberOfParticles, stochasticWeighting, stochasticDisplacement, tSeed, conditions.deltaTime);

		if( gDebug == 1 && stochasticWeighting != NULL)
        {
        //    conditions.currentTime = conditions.endTime+1;
            FILE *stochasticOutput = fopen("../bin/stochastic_matrix_output.txt","a");

            for(int i = 0; i < 6 * numberOfParticles; i++)
            {
                for(int j = 0; j < 6 * numberOfParticles; j++)
                {
                    fprintf(stochasticOutput, "%e\t", stochasticWeighting[i * 6 * numberOfParticles + j]);
                }
                fprintf(stochasticOutput, "\n");

            }
            fclose (stochasticOutput);
		}

        //
        // Include additional forces
        //

        force_torque_summation(additionalForces, generalisedCoordinates, 6 * numberOfParticles, forceList, numberOfForces, conditions);

        //
        // Calculate time step.
        //

        moving_on_routine(numberOfParticles, &conditions, diffusionMatrix, additionalForces, stochasticDisplacement, generalisedCoordinates, NULL);
        fprintf(output, "%lf\t", conditions.currentTime);
        for(int i = 0; i < 3 * numberOfParticles; i++)
        {
            fprintf(output, "%e\t", generalisedCoordinates[i]);
        }
        fprintf(output, "\n");

        conditions.currentTime+=conditions.deltaTime; // time step
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
