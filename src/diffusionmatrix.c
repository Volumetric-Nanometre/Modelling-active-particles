/*********************
* Date of creation 10/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
*
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <omp.h>

#include "diffusionmatrix.h"
#include "maths_functions.h"
#include "particles.h"

//
// Physucal constants
//

extern double gBoltzmannConst;
extern double gPi;

//
// Create the diffusion matrix
//

void diffusion_matrix_creation(int numberOfParticles, double *diffusionMatrix, double *stochasticWeighting, particleVariables *particles, environmentVariables *conditions)
{
    double tempTransMatrix[9] = {0} ;
    double tempRotatMatrix[9] = {0} ;
    double tempCouplMatrix[9] = {0} ;

    //
    // Scan through the particles and calculate the individual matrices
    //
    #pragma omp parallel for private(tempTransMatrix,tempRotatMatrix,tempCouplMatrix) collapse (2)
    for( int particleRow = 0; particleRow < numberOfParticles; particleRow++)
    {
        for( int particleColumn = 0 ; particleColumn < numberOfParticles; particleColumn++)
        {

            //
            // Create translational submatrix
            //
            translational_tensor_creation(tempTransMatrix, particles, conditions , particleRow, particleColumn );
            //
            // Create rotational submatrix
            //
            rotational_tensor_creation(tempRotatMatrix, particles, conditions, particleRow, particleColumn );
            //
            //Create coupled submatrix
            //
            translation_rotation_coupling_tensor_creation(tempCouplMatrix, particles, conditions, particleRow, particleColumn );
            //
            // Transfer to main diffusion matrix
            //
            for(int n = 0; n < 3; n ++)
            {
                for(int m = 0; m < 3; m ++)
                {
                    //
                    // Calculate the positions of the submatrices within the grand matrix
                    //
                    int translationPosition = (particleRow * 3 + n) * 6 * numberOfParticles + particleColumn * 3 + m;
                    int rotationPosition = ( (particleRow + numberOfParticles )* 3 + n) * 6 * numberOfParticles + (particleColumn + numberOfParticles) * 3 + m ;
                    int couplingPositionBottomLeft = ( (particleRow + numberOfParticles )* 3 + n) * 6 * numberOfParticles + particleColumn * 3 + m;
                    int couplingPositionTopRight = (particleRow* 3 + n) * 6 * numberOfParticles + (particleColumn + numberOfParticles) * 3 + m;
                    //
                    // Note the couplingPositionBottomLeft requires that the values be the negative of
                    // the couplingPositionTopRight values.
                    //
                    stochasticWeighting[translationPosition] = diffusionMatrix[ translationPosition ] = tempTransMatrix[n * 3 + m];
                    stochasticWeighting[rotationPosition] = diffusionMatrix[ rotationPosition ] = tempRotatMatrix[n * 3 + m];
                    stochasticWeighting[couplingPositionTopRight] = diffusionMatrix[ couplingPositionTopRight ] = tempCouplMatrix[n * 3 + m];
                    stochasticWeighting[couplingPositionBottomLeft] = diffusionMatrix[ couplingPositionBottomLeft ] = -tempCouplMatrix[n * 3 + m];
                }
            }
        }
    }

}


//
// Create 3 x 3 submatrices using the Oseen tensor rules
//
void translational_tensor_creation(double *tempMatrix, particleVariables *particles, environmentVariables *conditions, int i, int j)
{
    double stokesConstantProduct = ( gBoltzmannConst * conditions->temperature ) / ( gPi * conditions->viscosity);
    if(i==j)
    {
        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = kronecker_delta(n,m) * stokesConstantProduct / (6 * conditions->radius); // Over 6 for the self interaction terms


            }
        }
    }
    else
    {
        //
        // Calculate dyadic matrix elements by first calculating the dimensional
        // vectors, i.e x_ij = x_j - x_i
        //
        double dimensionalVector[3] = {0};
        dimensionalVector[0] = particles[j].x - particles[i].x;
        dimensionalVector[1] = particles[j].y - particles[i].y;
        dimensionalVector[2] = particles[j].z - particles[i].z;

        double absDistance = sqrt( pow(dimensionalVector[0],2) + pow(dimensionalVector[1],2) + pow(dimensionalVector[2],2) );

        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = ( kronecker_delta(n,m)+ (dimensionalVector[n] * dimensionalVector [m])/pow(absDistance,2) )
                                        *( stokesConstantProduct / 8 * absDistance);      //  Over 8 for interparticle interaction terms
            }

        }
    }

}

//
// Create the rotational diffusion 3x3 matrices
//

void rotational_tensor_creation(double *tempMatrix, particleVariables *particles, environmentVariables *conditions, int i, int j)
{
    double stokesConstantProduct = ( gBoltzmannConst * conditions->temperature ) / ( gPi * conditions->viscosity);
    if(i==j)
    {
        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = kronecker_delta(n,m) * stokesConstantProduct / (8 * pow(conditions->radius, 3) ); // Over 8 for the self interaction terms


            }
        }
    }
    else
    {
        //
        // Calculate dyadic matrix elements by first calculating the dimensional
        // vectors, i.e x_ij = x_j - x_i
        //
        double dimensionalVector[3] = {0};
        dimensionalVector[0] = particles[j].x - particles[i].x;
        dimensionalVector[1] = particles[j].y - particles[i].y;
        dimensionalVector[2] = particles[j].z - particles[i].z;

        double absDistance = sqrt( pow(dimensionalVector[0],2) + pow(dimensionalVector[1],2) + pow(dimensionalVector[2],2) );

        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] =( (dimensionalVector[n] * dimensionalVector [m]) / pow(absDistance,2) - kronecker_delta(n, m) )
                                        *( 3 * stokesConstantProduct / 16 * pow( absDistance, 3) );       //  Over 16 for interparticle interaction terms
            }

        }
    }
}


//
// Create the rotational translational coupled matrices
// Note T-R is the bottom left corner of the grand matrix and is the negative
// of T-R
//

void translation_rotation_coupling_tensor_creation(double *tempMatrix, particleVariables *particles, environmentVariables *conditions, int i, int j)
{
    double stokesConstantProduct = ( gBoltzmannConst * conditions->temperature ) / ( gPi * conditions->viscosity);
    if(i==j)
    {
        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = 0;
            }
        }
    }
    else
    {
        //
        // Calculate dyadic matrix elements by first calculating the dimensional
        // vectors, i.e x_ij = x_j - x_i
        //
        double dimensionalVector[3] = {0};
        dimensionalVector[0] = particles[j].x - particles[i].x;
        dimensionalVector[1] = particles[j].y - particles[i].y;
        dimensionalVector[2] = particles[j].z - particles[i].z;

        double absDistance = sqrt( pow(dimensionalVector[0],2) + pow(dimensionalVector[1],2) + pow(dimensionalVector[2],2) );

        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = ( dimensionalVector[n] / absDistance )*levi_civita_density(n,m)
                                    *( stokesConstantProduct / 8 * pow( absDistance, 2) );       //  Over 8 for interparticle interaction terms
            }

        }
    }
}
