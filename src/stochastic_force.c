
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
#include <gsl/gsl_linalg.h>
#include "stochastic_force.h"
#include "maths_functions.h"

#include <gsl/gsl_randist.h>

void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, gsl_rng *tSeed, double timestep)
{
	int N = 6 * numberOfParticles; // array is a linearized (6*N) by (6*N) array
	
	/*double value;
	printf("Initial array:\n");
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			value = stochasticWeighting[i*N + j];
			
			printf("%+1.1e", value);
			if (j<N-1) printf(" ");
		}
		printf("\n");
	}
	
	gsl_matrix *test = gsl_matrix_alloc(N, N);
	//gsl_matrix_set_identity(test);
	
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			gsl_matrix_set(test,i,j, stochasticWeighting[i*N + j]);
	
	int err = gsl_linalg_cholesky_decomp1(test);
	if (err != 0)
	{
		printf("Error occurred during cholesky decomposition\n");
		return;
	}
	else printf("Completed cholesky decomposition\n");*/
	

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

	for (int k = 0; k < N; k++)
	{
		for (int j = k + 1; j < N; j++)
		{
			stochasticWeighting[k*N+j] = 0; // iterates over upper triangle and sets all values to zero
		}
	}

	for (int i=0; i<N; i++)
	{
		stochasticDisplacement[i] = 0;
	}
	
	//gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

	double ran_num;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			//ran_num = gaussdev(&tSeed2);
			ran_num = gsl_ran_gaussian(tSeed, 1);
			printf("ran_num: %+1.5e\n", ran_num);
			stochasticDisplacement[i] += stochasticWeighting[i*N + j] * ran_num * sqrt(2*timestep);
		}
		printf("%+1.5e\n\n", stochasticDisplacement[i]);
	}
	printf("\nNEXT\n\n");
	
	//gsl_rng_free(r);
	
	/*gsl_matrix *array_gsl = gsl_matrix_alloc(N, N);
	
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
		{
			gsl_matrix_set(array_gsl, i, j, stochasticWeighting[i*N + j]);
			if (j>i)
			gsl_matrix_set(test, i, j, 0);
		}
			
	
	if (gsl_matrix_equal(test, array_gsl) == 1) printf("Cholesky decompositions agree!\n");
	else
	{
		printf("Cholesky decompositions do not agree\n");
		
		printf("Decomposed array:\n");
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
			{
				value = stochasticWeighting[i*N + j];
				if (value != gsl_matrix_get(array_gsl, i, j))
				{
					printf("Array and its GSL conversion don't match!\n");
					return;
				}
				
				printf("%+1.1e", value);
				if (j<N-1) printf(" ");
			}
			printf("\n");
		}
		
		printf("GSL decomposed array:\n");
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
			{
				printf("%+1.1e", gsl_matrix_get(test,i,j));
				if (j<N-1) printf(" ");
			}
			printf("\n");
		}
	}*/
}
