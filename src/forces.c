/*********************
* Date of creation 19/02/2018
* Authors: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* History
*
**************************************/

#include <stdlib.h>
#include <omp.h>
#include <math.h>

#include "forces.h"
#include "particles.h"



static void force_gravity(double *additionalForces, double numberOfCells, double mass);

static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, double numberOfCells, double radius);

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, double numberOfCells);

enum forces_available
{
    NONE ,
    GRAVITY ,
    VAN_DER_WAALS ,
    EXP_REPULSION
};

extern double gGrav;
//
// Tuning varibles for the exponential
// F = Aexp(-r/lamda)
// r= r2 -r1
//
static double gExpTuningA = -1;
static double gExpTuningLamda = 1;
//
// Calculate the sum of the forces and torques from the forces and torques wanted
//
void force_torque_summation(double *additionalForces,double *generalisedCoordinates, int numberOfCells, int *forceList, int numberOfForces,	environmentVariables conditions)
{
    //
    // Set entire array to zero
    //
    #pragma omp parallel for
    for(int i = 0; i < numberOfCells; i++)
    {
        additionalForces[i] = 0;
    }


    for(int i = 0; i < numberOfForces; i++)
    {
        switch(forceList[i])
        {
            case NONE : break;
            case GRAVITY : force_gravity(additionalForces, numberOfCells, conditions.mass);
            case VAN_DER_WAALS : force_van_der_waals(additionalForces, generalisedCoordinates, numberOfCells, conditions.radius);
            case EXP_REPULSION : force_exp_repulsion(additionalForces, generalisedCoordinates, numberOfCells);
            default : break;
        }
    }
}

static void force_gravity(double *additionalForces, double numberOfCells, double mass)
{
    //
    // Half the number of cells to ignore torques.
    // Step 3 each time to only hit the z axis
    // Add the forces onto the preexisting values
    //
    for(int i = 0; i < numberOfCells/2; i+=3)
    	additionalForces[i] -= mass * gGrav;
}

static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, double numberOfCells, double radius)
{
    //
    // get the number of particles.
    //
    int numberOfParticles = numberOfCells/6;
    // F = -32ArR^6/|r|^4(|r|^2 - 4R^2)^2
    // r = r2-r1
    // Add the forces onto the preexisting values
    //
    double x,y,z,r,forceConst;
    double hamakerCoeff=1; // This is A in the above eqn
    for(int i = 0; i < numberOfParticles; i++)
    {
        for(int j = 0; j < numberOfParticles; j++)
        {
            if(i == j)
            continue ;
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];

            r = sqrt(x * x + y * y + z * z);
            forceConst =  -32*hamakerCoeff*pow(radius,6) / ( 3 * pow(r,4)*pow( ( r*r - 4*radius*radius) , 2) );
            additionalForces[i * 3] += forceConst * x;
            additionalForces[i * 3 + 1] += forceConst * y;
            additionalForces[i * 3 + 2] += forceConst * z;
        }
    }
}

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, double numberOfCells)
{
    //
    // get the number of particles.
    //
    int numberOfParticles = numberOfCells/6;
    // F = Aexp(-r/lamda)
    // r = r2-r1
    // Add the forces onto the preexisting values
    //
    double x,y,z,r,forceConst;
    for(int i = 0; i < numberOfParticles; i++)
    {
        for(int j = 0; j < numberOfParticles; j++)
        {
            if(i == j)
            continue ;
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];

            r = sqrt(x * x + y * y + z * z);
            forceConst =  (gExpTuningA * exp(-r / gExpTuningLamda) ) / r;
            additionalForces[i * 3] += forceConst * x;
            additionalForces[i * 3 + 1] += forceConst * y;
            additionalForces[i * 3 + 2] += forceConst * z;
        }
    }

}
