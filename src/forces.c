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



static void force_gravity(double *additionalForces, int numberOfCells, double mass);

static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double radius);

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, int numberOfCells);

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells);

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField);

extern int gDebug;

extern double gGrav;
extern double gPi;
//
// Tuning varibles for the exponential
// F = Aexp(-r/lamda)
// r= r2 -r1
//
static double gExpTuningA = 7.176291667E-11 * 6E11;
static double gExpTuningLamda = 4E5;
//
// Calculate the sum of the forces and torques from the forces and torques wanted
//
void force_torque_summation(double *additionalForces,double *generalisedCoordinates, int numberOfCells, int *forceList, int numberOfForces,	environmentVariables conditions, field_t drivingField)
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
            case GRAVITY : force_gravity(additionalForces, numberOfCells, conditions.mass); break;
            case VAN_DER_WAALS : force_van_der_waals(additionalForces, generalisedCoordinates, numberOfCells, conditions.radius); break;
            case EXP_REPULSION : force_exp_repulsion(additionalForces, generalisedCoordinates, numberOfCells); break;
			case ALIGN_TORQUE : alignment_torque(additionalForces, generalisedCoordinates, numberOfCells); break;
			case DRIVING_FIELD : driving_force(additionalForces, generalisedCoordinates, numberOfCells, drivingField); break;
            default : break;
        }
    }
}

static void force_gravity(double *additionalForces, int numberOfCells, double mass)
{
    //
    // Half the number of cells to ignore torques.
    // Step 3 each time to only hit the z axis
    // Add the forces onto the preexisting values
    //
    for(int i = 0; i < numberOfCells/2; i+=3)
    	additionalForces[i] -= mass * gGrav;
}

static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double radius)
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
    double hamakerCoeff=2.5E-19; // This is A in the above eqn, constant for unretarded gold
    for(int i = 0; i < numberOfParticles; i++)
    {
        for(int j = 0; j < numberOfParticles; j++)
        {
            //
            // Ignore self interaction
            //
            if(i == j)
            continue ;
            //
            // Calculate the vector components of r
            //
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];
            //
            // Calculate |r|
            //
            r = sqrt(x * x + y * y + z * z);
            //
            // Calculate  -32AR^6 / 3|r|^4(|r|^2 - 4R^2)^2
            //
            forceConst =  -32*hamakerCoeff*pow(radius,6) / ( 3 * pow(r,4)*pow( ( r*r - 4*radius*radius) , 2) );
            //
            // Vectorise
            //
            additionalForces[i * 3] += forceConst * x;
            additionalForces[i * 3 + 1] += forceConst * y;
            additionalForces[i * 3 + 2] += forceConst * z;
        }
    }
}

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, int numberOfCells)
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
            //
            // Ignore self interaction
            //
            if(i == j)
            continue ;
            //
            // Calculate the vector components of r
            //
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];
            //
            // Calculate |r|
            //
            r = sqrt(x * x + y * y + z * z);
            //
            // Calculate 1/R * Aexp(-|r|/lamda)
            //
            forceConst =  gExpTuningA * exp(-gExpTuningLamda * r) / r;
            //
            // Vectorise
            //
            additionalForces[i * 3] += forceConst * x;
            additionalForces[i * 3 + 1] += forceConst * y;
            additionalForces[i * 3 + 2] += forceConst * z;
        }
    }

}

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells)
{
	double totalAlpha = 0;
	double totalBeta = 0;
	
	int torqueOffset = 3*numberOfCells;
	
	double meanAlpha, meanBeta;
	
	double forceConst = 1;
	
	// Calculate average angles
	for (int i=0; i<numberOfCells/2; i++)
	{
		totalAlpha += generalisedCoordinates[torqueOffset + 3*i + 0];
		totalBeta += generalisedCoordinates[torqueOffset + 3*i + 1];
	}
	
	meanAlpha = totalAlpha/numberOfCells;
	meanBeta = totalBeta/numberOfCells;
	
	// Calculate torques on each particle according to their alignment with average angle
	for (int i=0; i<numberOfCells/6; i++)
	{
		additionalForces[torqueOffset + 3*i + 0] += forceConst * sin(meanAlpha - generalisedCoordinates[torqueOffset + 3*i + 0]);
		additionalForces[torqueOffset + 3*i + 1] += forceConst * sin(meanBeta - generalisedCoordinates[torqueOffset + 3*i + 1]);
	}
	
}

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField)
{
	
	for (int i=0; i<numberOfCells/6; i++)
	{
		additionalForces[3*i + 0] +=  drivingField.mag * cos(drivingField.alpha - generalisedCoordinates[numberOfCells/2 + 3*i + 0]);
		additionalForces[3*i + 1] +=  drivingField.mag * cos(drivingField.beta - generalisedCoordinates[numberOfCells/2 + 3*i + 1]);
	}
	
}
