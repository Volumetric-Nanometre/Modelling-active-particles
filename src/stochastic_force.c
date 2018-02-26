
/*********************
* Date of creation 17/10/2017
* Author: Oliver Hinds
* Contact:
**************************************
* History
**************************************/

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include "stochastic_force.h"
#include "maths_functions.h"


extern int gDebug;


void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, gsl_rng *tSeed, double timestep){
	
	int N = 6 * numberOfParticles; // array is a linearized (6*N) by (6*N) array
	
	for (int k = 0; k < N; k++) // iterates over diagonals
	{
		stochasticWeighting[k*(N+1)] = sqrt(stochasticWeighting[k*(N+1)]); // square roots diagonals

		for (int i = k + 1; i < N; i++) // iterates over the elements in the column below the diagonal
		{
			stochasticWeighting[i*N + k] = stochasticWeighting[i*N + k] / stochasticWeighting[k*(N+1)];
		}


		// iterate over lower triangle subtended by (k,k) element
		for (int j = k + 1; j < N; j++) // iterates over the columns j>k
		{
			for (int i = j; i < N; i++) // iterates over the rows i>k
			{
				stochasticWeighting[i*N + j] = stochasticWeighting[i*N + j] - (stochasticWeighting[i*N + k] * stochasticWeighting[j*N + k]);
			}
		}
    }
	
	// iterate over upper triangle and set all values to zero
	for (int k = 0; k < N; k++)
	{
		for (int j = k + 1; j < N; j++)
		{
			stochasticWeighting[k*N+j] = 0;
		}
	}

	double ran_num;
	for (int i = 0; i < N; i++)
	{
		stochasticDisplacement[i] = 0;
		for (int j = 0; j < N; j++)
		{
			ran_num = gsl_ran_gaussian(tSeed, 1);
			stochasticDisplacement[i] += stochasticWeighting[i*N + j] * ran_num * sqrt(2*timestep);
		}
		//if (gDebug == 1) printf("%+1.5e\n\n", stochasticDisplacement[i]);
	}
	//if (gDebug == 1) printf("\nNEXT\n\n");
}
