
/*********************
* Date of creation 17/10/2017
* Author: Oliver Hinds
* Contact:
* Other Authors: N/A
**************************************
* History
**************************************/

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "stochastic_force.h"
#include "maths_functions.h"

void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, long int tSeed)
{
	int N = 6 * numberOfParticles; // array is a linearized (6*N) by (6*N) array


	int M = 6*numberOfParticles;
   gsl_matrix  *A=gsl_matrix_alloc(M,M);

   for(int i =0; i<M; i++)
   {
       for(int j = 0; j< M; j++)
       {
           gsl_matrix_set( A,i,  j, stochasticWeighting[i*M+j]);
       }
   }

   gsl_linalg_cholesky_decomp1( A);
   getchar();

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



	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			stochasticDisplacement[i] += stochasticWeighting[i*N + j] * guassdev(&tSeed); // X[j] is some random gaussian

		}
	}

}
