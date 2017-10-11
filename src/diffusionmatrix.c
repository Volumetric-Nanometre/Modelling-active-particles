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

double *diffusion_matrix_creation(int numberOfParticles, particleVariables *particles, double temperature, double viscosity, double radius)
{

    double *diffusionMatrix = NULL ;
    double tempMatrix[9] = {0} ;

    //
    // Allocate and check memory
    //
    diffusionMatrix = calloc( pow(6 * numberOfParticles, 2 ), sizeof ( *diffusionMatrix ));

    if( diffusionMatrix == NULL )
    {
        printf("-Error %d : %s\n", errno, strerror( errno ) );
        return NULL;
    }

    //
    // Scan through the particles and calculate the individual Oseen matrices
    //
    for( int particleRow = 0; particleRow < 3*numberOfParticles; particleRow+=3)
    {
        for( int particleColumn = 0 ; particleColumn < 3*numberOfParticles; particleColumn+=3)
        {
            //
            // Create Oseen matrix
            //
            translational_tensor_creation(tempMatrix, particles, temperature, viscosity, radius, particleRow, particleColumn );

            //
            // Transfer to main diffusion matrix
            //
            for(int n = 0; n < 3; n ++)
            {
                for(int m = 0; m < 3; m ++)
                {
                    diffusionMatrix[ (particleRow + n)*3*numberOfParticles +particleColumn+m ] = tempMatrix[n * 3 + m];
                }
            }

        }
    }

    for( int particleRow = 3*numberOfParticles; particleRow < 6*numberOfParticles; particleRow+=3)
    {
        for( int particleColumn = 3*numberOfParticles ; particleColumn < 6*numberOfParticles; particleColumn+=3)
        {
            //
            // Create Oseen matrix
            //
            rotational_tensor_creation(tempMatrix, particles, temperature, viscosity, radius, particleRow, particleColumn );

            //
            // Transfer to main diffusion matrix
            //
            for(int n = 0; n < 3; n ++)
            {
                for(int m = 0; m < 3; m ++)
                {
                    diffusionMatrix[ (particleRow + n)*3*numberOfParticles +particleColumn+m ] = tempMatrix[n * 3 + m];
                }
            }

        }
    }

    return diffusionMatrix;
}


//
// Create 3 x 3 submatrices using the Oseen tensor rules
//
void translational_tensor_creation(double *tempMatrix, particleVariables *particles, double temperature, double viscosity, double radius, int i, int j)
{
    double stokesConstantProduct = ( gBoltzmannConst * temperature ) / ( gPi * viscosity);
    if(i==j)
    {
        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = kronecker_delta(n,m) * stokesConstantProduct / (6* radius); // Over 6 for the self interaction terms


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

void rotational_tensor_creation(double *tempMatrix, particleVariables *particles, double temperature, double viscosity, double radius, int i, int j)
{
    double stokesConstantProduct = ( gBoltzmannConst * temperature ) / ( gPi * viscosity);
    if(i==j)
    {
        for(int n = 0; n < 3; n ++)
        {
            for(int m = 0; m < 3; m ++)
            {
                tempMatrix[n * 3 + m] = kronecker_delta(n,m) * stokesConstantProduct / (8 * pow(radius, 3) ); // Over 8 for the self interaction terms


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

void translation_rotation_coupling_tensor_creation(double *tempMatrix, particleVariables *particles, double temperature, double viscosity, double radius, int i, int j)
{
    double stokesConstantProduct = ( gBoltzmannConst * temperature ) / ( gPi * viscosity);
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
