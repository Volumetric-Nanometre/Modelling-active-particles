
/*********************
* Date of creation 09/10/2017
* Author: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* History
*
* 09/10/2017 -Created particle read in function to allocate memory and
*             bring in inputs.
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "initial_finalisation.h"
#include "particles.h"
#include "maths_functions.h"

extern double gPi;

int particle_read_in(particleVariables **particles)
{
    int numberOfParticles=0;

    particleVariables* initParticles=NULL;

    //
    // Open and check input file
    //
    FILE* particleInput=fopen("../bin/particleInput.txt","r");

    if(particleInput == NULL)
    {
        printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
        numberOfParticles=-errno;
    }
    else
    {

        //
        // Get number of particles. If Illegal value, close file and return -1
        //

        if(fscanf(particleInput,"%d\n",&numberOfParticles)!=1 || numberOfParticles<=0 )
        {
            printf("- Error : Number of particles is invalid \n");
            numberOfParticles=0;
            fclose(particleInput);
        }
        else
        {
            //
            // Allocate memory for input then read until
            // # of particles reached or the last whole line is read.
            //

            if((initParticles=calloc(numberOfParticles,sizeof(*initParticles))) == NULL)
            {
                printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
                fclose (particleInput);
                numberOfParticles=-errno;
            }

            int i=0;

            while(i<numberOfParticles && fscanf(particleInput,"%lf%lf%lf%lf%lf%lf\n",&initParticles[i].x,&initParticles[i].y,&initParticles[i].z,
                                                                                    &initParticles[i].alpha,&initParticles[i].beta,&initParticles[i].gamma)==6)
            {
                i++;
            }

        }
    }
    fclose(particleInput);
    *particles=initParticles;
    return numberOfParticles;
}

int generate_particle_data(int numberOfParticles, particleVariables **particles, gsl_rng *tSeed, double xMax, double yMax, double zMax)
{
	particleVariables *initParticles = calloc(numberOfParticles, sizeof(*initParticles));

	for (int i=0; i<numberOfParticles; i++)
	{
		initParticles[i].x = xMax * gsl_rng_uniform(tSeed);
		initParticles[i].y = yMax * gsl_rng_uniform(tSeed);
		initParticles[i].z = zMax * gsl_rng_uniform(tSeed);
		initParticles[i].alpha = 2*gPi * gsl_rng_uniform(tSeed);
		initParticles[i].beta = 2*gPi * gsl_rng_uniform(tSeed);
		initParticles[i].gamma = 2*gPi * gsl_rng_uniform(tSeed);
	}

	*particles = initParticles;

	return 0;
}

//
// Takes the spatial and rotational coordinates of the particles and creates
// a generalised coordinate to represent every particle.
//

double *generalised_coordinate_creation(int numberOfParticles, particleVariables *particles)
{
    double *generalisedCoordinates = NULL ;
    //
    // Allocate and check memory
    //
    generalisedCoordinates = calloc( 6 * numberOfParticles, sizeof  *generalisedCoordinates );

    if( generalisedCoordinates == NULL )
    {
        printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
        return NULL;
    }

    //
    // Allocate the coordinates
    //

    for(int i = 0; i < numberOfParticles; i ++)
    {
        generalisedCoordinates[ i * 3 ] = particles[i].x;
        generalisedCoordinates[ i * 3 + 1 ] = particles[i].y;
        generalisedCoordinates[ i * 3 + 2 ] = particles[i].z;
        generalisedCoordinates[ ( i + numberOfParticles ) * 3 ] = particles[i].alpha;
        generalisedCoordinates[ ( i + numberOfParticles ) * 3 + 1 ] = particles[i].beta;
        generalisedCoordinates[ ( i + numberOfParticles ) * 3 + 2 ] = particles[i].gamma;
    }

    return generalisedCoordinates;
}
