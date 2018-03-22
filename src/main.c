/*********************
* Date of creation 09/10/2017
* Authors: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* Change History
**************************************/
//
// GNU and external libraries
//
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <omp.h>
#include "mpi.h"
//
// Bespoke libraries
//
#include "particles.h"
#include "diffusionmatrix.h"
#include "stochastic_force.h"
#include "moving_on.h"
#include "forces.h"
#include "initial_finalisation.h"

#define MASTER 0
#define TOP_SLAVE 1


double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;
double gGrav = 9.80665; // m s^-2

int gDebug = 0;
int gSerial = 0;
int gNumOfthreads;
int gNumOfNodes = 1;

enum MESSAGE_TAGS
{
    NUM_OF_PARTICLES,
    NUM_OF_FORCES ,         // Number of forces to expect
    FORCE_LIST,             // Which forces to use
    COORDINATES ,           // Corordinates to use
    ADDITIONAL_FORCES       // Total forces
};

int main(int argc, char *argv[])
{
    //
    // MPI Initilisation
    //
    int taskid;
    int MPI_error = 0;
    MPI_Status status;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &gNumOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);


	//
	// Create conditions variable
	//
	environmentVariables conditions;
	//
	// Initilise conditions with boilerplate varaibles
	//
	boilerplate_variables(&conditions);
	//
	// Read in cmd line arguments and adjust conditions as neccessary
	//
	cmd_line_read_in(argc, argv, &conditions);
    //
    // Master process environment
    //
    // Create driving field
    field_t drivingField;
    drivingField.mag = 1E-10;
    drivingField.alpha = 0;
    drivingField.beta = gPi/2;


    //
    //  Choose forces to be included
    //
    int numberOfForces = 2; // must be at least 1, with the force none chosen
    //
    // Copy of enum to understand force forceList
    //
    //enum forces_available
    //{
    //    NONE ,
    //    GRAVITY ,
    //    VAN_DER_WAALS ,
    //    EXP_REPULSION ,
    //	  ALIGN_TORQUE ,
    //	  DRIVING_FIELD,
    //	  POLAR_DRIVING_FORCE
    //};


    int forceList[4] = {VAN_DER_WAALS,EXP_REPULSION, POLAR_DRIVING_FORCE, ALIGN_TORQUE};

    if( MASTER == taskid )
    {
        //
    	//Open files for output
    	//
        FILE *output = fopen("../bin/output.csv","w");
        FILE *angle_output = fopen("../bin/angle_output.csv","w");
    	FILE *forces_output = fopen("../bin/forces_output.csv","w");
        if(output == NULL || angle_output == NULL || forces_output == NULL)
        {
            printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
            return -errno;
        }
    	//
    	// Initilise random variables
    	//
    	gsl_rng **rndarray=rand_array_allocation();

    	if(rndarray == NULL)
    	{
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
    		return -1;
    	}

        //
        // Initilise generalised coordinates
        //

    	double *generalisedCoordinates = generalised_coordinate_initilisation(&conditions,rndarray);

    	if(generalisedCoordinates == NULL)
    	{
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
    		return -1;
    	}

        int vectorSize = 6 * conditions.numberOfParticles;
        //
        // Send number of particles to slaves
        //
        if(gNumOfNodes > 1)
        {
            MPI_Send(&conditions.numberOfParticles, 1, MPI_INT, TOP_SLAVE, NUM_OF_PARTICLES, MPI_COMM_WORLD);
        }

        //---------------------------- DEBUG ------------------------------//
        //
        // Prints the generalisedCoordinates to a file for inspection
        //
        if( gDebug == 1 && generalisedCoordinates != NULL)
        {
            FILE *genCoordOutput = fopen("../bin/genCoord_output.txt","w");

            for(int i = 0; i < vectorSize; i++)
            {
                fprintf(genCoordOutput, "%e\n", generalisedCoordinates[i]);
            }
            fclose (genCoordOutput);
        }
        //--------------------------- END ---------------------------------//


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

        diffusionMatrix = calloc( pow(vectorSize, 2), sizeof *diffusionMatrix) ;
        stochasticWeighting = calloc( pow( vectorSize, 2), sizeof *stochasticWeighting);
        stochasticDisplacement = calloc( vectorSize, sizeof *stochasticDisplacement);
        additionalForces = calloc( vectorSize, sizeof *additionalForces);

        if(  diffusionMatrix==NULL  || stochasticWeighting==NULL || stochasticDisplacement==NULL || additionalForces==NULL)
        {
    		free_memory(6,diffusionMatrix, generalisedCoordinates, stochasticWeighting, stochasticDisplacement,additionalForces);
    		diffusionMatrix = generalisedCoordinates = stochasticWeighting = stochasticDisplacement = additionalForces = NULL ;
            printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
            return -errno;
        }


        //
        // Loop through time, output each time step to a file.
        //
        int loop = 0;
        int count = 0;
        int maxLoop = conditions.endTime/(double)conditions.deltaTime;

        double progTime = omp_get_wtime();

        while(conditions.currentTime<=conditions.endTime)
        {
            //
            // Create diffusion matrix
            //
            if(gNumOfNodes > 1)
            {
                MPI_Send(&generalisedCoordinates[0], vectorSize, MPI_DOUBLE, TOP_SLAVE, COORDINATES, MPI_COMM_WORLD);
            }

    	    diffusion_matrix_creation( conditions.numberOfParticles, diffusionMatrix, stochasticWeighting, generalisedCoordinates, &conditions);

    	    //---------------------------- DEBUG------------------------------//
    	    //
    	    // Prints the diffusionMatrix to a file for inspection
    	    //
    	    if( gDebug == 1 && diffusionMatrix != NULL)
    	    {
    	        //conditions.currentTime = conditions.endTime+1;
    	        FILE *matrixOutput = fopen("../bin/matrix_output.txt","w");

    	        for(int i = 0; i < 6 * conditions.numberOfParticles; i++)
    	        {
    	            for(int j = 0; j < 6 * conditions.numberOfParticles; j++)
    	            {
    	                fprintf(matrixOutput, "%e\t", diffusionMatrix[i * 6 * conditions.numberOfParticles + j]);
    	            }
    	            fprintf(matrixOutput, "\n");

    	        }
    	        fclose (matrixOutput);
    	    }
    	    //---------------------------END---------------------------------//


    	    //
    	    // Create the stochastic displacement
    	    //

    	    stochastic_displacement_creation( conditions.numberOfParticles, stochasticWeighting, stochasticDisplacement, rndarray, conditions.deltaTime);

    		if( gDebug == 1 && stochasticWeighting != NULL)
    	    {
    	    	//conditions.currentTime = conditions.endTime+1;
    	        FILE *stochasticOutput = fopen("../bin/stochastic_matrix_output.txt","w");

    	        for(int i = 0; i < 6 * conditions.numberOfParticles; i++)
    	        {
    	            for(int j = 0; j < 6 * conditions.numberOfParticles; j++)
    	            {
    	                fprintf(stochasticOutput, "%e\t", stochasticWeighting[i * 6 * conditions.numberOfParticles + j]);
    	            }
    	            fprintf(stochasticOutput, "\n");
    	        }
    	        fclose (stochasticOutput);
    		}


    		//
    		// Include additional forces (only calculates if its the only node)
    		//
            if(gNumOfNodes <= 1)
            {
                force_torque_summation(additionalForces, generalisedCoordinates, 6 * conditions.numberOfParticles, forceList, numberOfForces, conditions, drivingField);
            }
            else
            {
                MPI_Recv(&additionalForces[0], vectorSize, MPI_DOUBLE, TOP_SLAVE, ADDITIONAL_FORCES, MPI_COMM_WORLD, &status);
            }

            //
            // Calculate time step.
            //
            moving_on_routine(conditions.numberOfParticles, &conditions, diffusionMatrix, additionalForces, stochasticDisplacement, generalisedCoordinates, NULL);
            if(loop%100 == 0)
            {
    			int angle_offset = 3*conditions.numberOfParticles;
                fprintf(output, "%e, ", conditions.currentTime);
                fprintf(angle_output, "%e, ", conditions.currentTime);
    			fprintf(forces_output, "%e,",conditions.currentTime);
                fflush(output);
                fflush(angle_output);
                fflush(forces_output);
                for(int i = 0; i < 3 * conditions.numberOfParticles; i++)
                {
                    fprintf(output, "%e", generalisedCoordinates[i]);
    				fprintf(forces_output, "%e", additionalForces[i]);
                    fprintf(angle_output, "%e", fmod(generalisedCoordinates[angle_offset + i],2*gPi));
                    fflush(output);
                    fflush(angle_output);
                    fflush(forces_output);
    				if (i < 3*conditions.numberOfParticles - 1)
                    {
    	                fprintf(output, ", ");
    	                fprintf(angle_output, ", ");
    					fprintf(forces_output, ", ");
                        fflush(output);
                        fflush(angle_output);
                        fflush(forces_output);
                    }
                }
                fprintf(output, "\n");
                fprintf(angle_output, "\n");
    			fprintf(forces_output, "\n");
                fflush(output);
                fflush(angle_output);
                fflush(forces_output);
            }

    		loop ++;
            conditions.currentTime+=conditions.deltaTime; // time step
    		if((maxLoop/10)*count == loop)
    		{
    			//printf("%d%%\n",count*10);
    			count++;
    		}
            loop ++;
        }
    	progTime = omp_get_wtime() - progTime;

    	printf("Run time %lf s\n",progTime);

        //
        // Free memory
        //

        fclose(output);
    	fclose(angle_output);
    	fclose(forces_output);

    	free_memory(5,diffusionMatrix, generalisedCoordinates, stochasticWeighting, stochasticDisplacement, additionalForces);
    	diffusionMatrix = generalisedCoordinates = stochasticWeighting = stochasticDisplacement = additionalForces = NULL ;
    }

    if( TOP_SLAVE == taskid )
    {
        //
        // Recieve number of particles and update conditions
        //
        MPI_Recv(&conditions.numberOfParticles, 1, MPI_INT, MASTER, NUM_OF_PARTICLES, MPI_COMM_WORLD, &status);

        int vectorSize = 6 * conditions.numberOfParticles;

        //
        // Allocate memory required for the program.
        // Requires: additional forces, generalised coordinates
        //
        double *additionalForces = NULL;
        double *generalisedCoordinates = NULL;

        generalisedCoordinates = calloc( vectorSize, sizeof *generalisedCoordinates);
        additionalForces = calloc( vectorSize, sizeof *additionalForces);

        if( generalisedCoordinates==NULL || additionalForces==NULL)
        {
            free_memory(2, generalisedCoordinates,additionalForces);
            generalisedCoordinates = additionalForces = NULL ;
            printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
            return -errno;
        }

        while(conditions.currentTime<=conditions.endTime)
        {
            //
            // Recieve coordinates
            //
            MPI_Recv(&generalisedCoordinates[0], vectorSize, MPI_DOUBLE, MASTER, COORDINATES, MPI_COMM_WORLD, &status);
    		//
    		// Calculate additional forces
    		//
            force_torque_summation(additionalForces, generalisedCoordinates, 6 * conditions.numberOfParticles, forceList, numberOfForces, conditions, drivingField);
            //
            // Send new forces
            //
            MPI_Send(&additionalForces[0], vectorSize, MPI_DOUBLE, MASTER, ADDITIONAL_FORCES, MPI_COMM_WORLD);

            conditions.currentTime+=conditions.deltaTime; // time step
    	}
        free_memory(2, generalisedCoordinates,additionalForces);
        generalisedCoordinates = additionalForces = NULL ;

    }
    MPI_Finalize();
    return 0;
}
