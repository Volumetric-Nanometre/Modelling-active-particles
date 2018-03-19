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
#include <gsl/gsl_rng.h>

#include "particles.h"
#include "diffusionmatrix.h"
#include "stochastic_force.h"
#include "moving_on.h"
#include "forces.h"

#include <gsl/gsl_rng.h>

double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;
double gGrav = 9.80665; // m s^-2

int gDebug = 0;
int gSerial = 0;

int main(int argc, char *argv[])
{
	int numberOfParticles = 0;
	
	double xMax = 1E-7;
	double yMax = 1E-7;
	double zMax = 1E-7;
	
	//
    // Allocate the environmental conditions and nano particle
    // characteristics
    //

	environmentVariables conditions;
    conditions.temperature = 298; // K
    conditions.viscosity = 8.9E-4; //N m^-2 s
    conditions.radius = 50E-9; // m
    conditions.currentTime = 0; // Seconds
    conditions.deltaTime = 1E-5; // Seconds
    conditions.endTime = 15; // Seconds
	conditions.mass = (4/3) * gPi * pow(conditions.radius,3)*19320; // kg - density of gold
	
	if (argc > 1)
	{
		for (int i=0; i<argc; i++)
		{
			if(strstr(argv[i],"-debug") != NULL) gDebug = 1;
			else if(strstr(argv[i],"-serial") != NULL)	gSerial = 1;
			else if (strstr(argv[i],"-num") != NULL || strstr(argv[i],"-n") != NULL)
			{
				if (sscanf(argv[i+1],"%d", &numberOfParticles) != 1)
				{
					printf("Invalid number of particles\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-cube") != NULL)
			{
				double temp_num;
				if (strstr(argv[i+1],"R") != NULL)
				{
					printf("R detected\n");
					if (sscanf(argv[i+1],"%lfR", &temp_num) != 2)
					{
						xMax = temp_num * conditions.radius;
						yMax = temp_num * conditions.radius;
						zMax = temp_num * conditions.radius;
					}
					else
					{
						printf("Invalid maximum dimension values\n");
						return -1;
					}
				}
				else if (sscanf(argv[i+1],"%lf", &temp_num) == 1)
				{
					xMax = temp_num;
					yMax = temp_num;
					zMax = temp_num;
				}
				else
				{
					printf("Invalid maximum dimension values\n");
					return -1;
				}
				printf("Set dimensions to %1.1em\n", xMax);
			}
			else if (strstr(argv[i],"-x") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &xMax) != 1)
				{
					printf("Invalid maximum x-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-y") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &yMax) != 1)
				{
					printf("Invalid maximum y-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-z") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &zMax) != 1)
				{
					printf("Invalid maximum z-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-t") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions.endTime) != 1)
				{
					printf("Invalid duration\n");
					return -1;
				}
				else printf("Simulation duration set to %1.1e seconds\n", conditions.endTime);
			}
			else if (strstr(argv[i],"-dt") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions.deltaTime) != 1)
				{
					printf("Invalid timestep\n");
					return -1;
				}
				else printf("Simulation timestep set to %1.1e seconds\n", conditions.deltaTime);
			}
		}
		if (gDebug == 1 && gSerial == 1) printf("Debug & serial modes active\n");
		else if (gDebug == 1 && gSerial == 0) printf("Debug mode active\n");
		else if (gDebug == 0 && gSerial == 1) printf("Serial mode active\n");

	}
	

    FILE *output = fopen("../bin/output.csv","w");
    FILE *angle_output = fopen("../bin/angle_output.csv","w");
	FILE *forces_output = fopen("../bin/forces_output.csv","w");
    if(output == NULL)
    {
        printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
        return -errno;
    }
	
	
	/*time_t tSeed1;
	time(&tSeed1);
	long int tSeed = -1*(long int) tSeed1;*/
	gsl_rng *tSeed = gsl_rng_alloc(gsl_rng_mt19937);
	
	particleVariables* particles = NULL;
	
	if (numberOfParticles == 0)
	{
	    //
	    // Call function to read in particle data
	    //
	    if( (numberOfParticles = particle_read_in(&particles)) <= 0)
	    {
	        getchar();
	        return numberOfParticles;
	    }
	
	    printf("Data read in success\n" );
	}
	else
	{
		generate_particle_data(numberOfParticles, &particles, tSeed, xMax, yMax, zMax);
		
		printf("Generated particle data\n");
	}
    
	
	// Create driving field
	field_t drivingField;
	drivingField.mag = 1E-10;
	drivingField.alpha = 0;
	drivingField.beta = gPi/2;

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
    //  Choose forces to be included
    //

    int numberOfForces = 2; // must be at least 1, with the force none chosen

	/* Forces available:
		NONE,
		GRAVITY,
		VAN_DER_WAALS,
		EXP_REPULSION,
		ALIGN_TORQUE,
		DRIVING_FIELD
	*/
    int forceList[2] = {VAN_DER_WAALS, EXP_REPULSION,};

    //
    // Loop through time, output each time step to a file.
    //
    int loop = 0;
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
            FILE *matrixOutput = fopen("../bin/matrix_output.txt","w");

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
        	//conditions.currentTime = conditions.endTime+1;
            FILE *stochasticOutput = fopen("../bin/stochastic_matrix_output.txt","w");

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

        force_torque_summation(additionalForces, generalisedCoordinates, 6 * numberOfParticles, forceList, numberOfForces, conditions, drivingField);

        //
        // Calculate time step.
        //

        moving_on_routine(numberOfParticles, &conditions, diffusionMatrix, additionalForces, stochasticDisplacement, generalisedCoordinates, NULL);
        //if(loop%10000 == 0)
        {
			int angle_offset = 3*numberOfParticles;
            fprintf(output, "%e, ", conditions.currentTime);
            fprintf(angle_output, "%e, ", conditions.currentTime);
			fprintf(forces_output, "%e,",conditions.currentTime);
            for(int i = 0; i < 3 * numberOfParticles; i++)
            {
                fprintf(output, "%e", generalisedCoordinates[i]);
				fprintf(forces_output, "%e", additionalForces[i]);
                fprintf(angle_output, "%e", fmod(generalisedCoordinates[angle_offset + i],2*gPi));
				if (i < 3*numberOfParticles - 1)
	                fprintf(output, ", ");
	                fprintf(angle_output, ", ");
					fprintf(forces_output, ", ");
            }
            fprintf(output, "\n");
            fprintf(angle_output, "\n");
			fprintf(forces_output, "\n");
        }
		
		loop ++;
        conditions.currentTime+=conditions.deltaTime; // time step
    }

    //
    // Free memory
    //

    fclose (output);
	fclose(angle_output);
	fclose(forces_output);

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
