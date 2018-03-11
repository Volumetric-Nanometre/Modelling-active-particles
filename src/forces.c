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

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells, environmentVariables conditions);

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField);

extern int gDebug;

extern double gGrav;
extern double gPi;
//
// Tuning varibles for the exponential
// F = Aexp(-r/lamda)
// r= r2 -r1
//
static double gExpTuningA = 6.29549E-12;
static double gExpTuningB = 3.519E7;


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
			case ALIGN_TORQUE : alignment_torque(additionalForces, generalisedCoordinates, numberOfCells, conditions); break;
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
//    #pragma omp parallel for
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
//    #pragma omp parallel for private(x,y,z,r,forceConst)
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
            additionalForces[j * 3] += forceConst * x;
            additionalForces[j * 3 + 1] += forceConst * y;
            additionalForces[j * 3 + 2] += forceConst * z;
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
    double expTuningA_B = gExpTuningA/gExpTuningB;
//    #pragma omp parallel for private(x,y,z,r,forceConst)
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
            forceConst =  gExpTuningA* exp(-gExpTuningB * r) / r;
            //
            // Vectorise
            //
            additionalForces[j * 3] += forceConst * x;
            additionalForces[j * 3 + 1] += forceConst * y;
            additionalForces[j * 3 + 2] += forceConst * z;
        }
    }

}

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells, environmentVariables conditions)
{
	int numberOfParticles = numberOfCells/6;
	int rotOffset = numberOfCells/2;
<<<<<<< HEAD
	
	// Calculated from rotational version of Langevin equation, substituting that the average angular displacement per timestep is pi/2
	// Torque T= pi*I/dt^2 for change in angle ~pi/2 per timestep
	// Moment of intertia I = (2/5)MR^2 for solid sphere of radius R and mass M
	double forceConst = gPi * 0.4 * conditions.mass * pow(conditions.radius, 2) * conditions.radius / pow(conditions.deltaTime, 2);
	/*															  |						  |
															Kept separate simply for clarity, might
														  be more comptuationally efficient to combine.
		 The second radius multiplier comes from using the inverse distance mutliplier in the force equation
		 	- the result is that the force is normalised based on particles separated by a distance of the same order as their radius
	*/
	
=======

	double forceConst = 1E-7;

>>>>>>> refs/remotes/origin/Michael
	double totalX = 0;
	double totalY = 0;
	double totalZ = 0;
	double totalAlpha = 0;
	double totalBeta = 0;

	double meanX, meanY, meanZ;
	double meanAlpha, meanBeta, difAlpha, difBeta;
<<<<<<< HEAD
	
	double dist,distMul;
	
	// Sum position and angles
=======

	double distMul;

	// Calculate average position and angles
>>>>>>> refs/remotes/origin/Michael
	for (int i=0; i<numberOfParticles; i++)
	{
		totalX += generalisedCoordinates[3*i + 0];
		totalY += generalisedCoordinates[3*i + 1];
		totalZ += generalisedCoordinates[3*i + 2];
<<<<<<< HEAD
		
		// Summing only alpha and beta angles
		totalAlpha += generalisedCoordinates[rotOffset + 3*i + 0];
		totalBeta += generalisedCoordinates[rotOffset + 3*i + 1];
	}
	
	// Calculate average positions
	meanX = totalX/numberOfParticles;
	meanY = totalY/numberOfParticles;
	meanZ = totalZ/numberOfParticles;
	
	// Calculate average angles
=======

		totalAlpha += generalisedCoordinates[rotOffset + 3*i + 0];
		totalBeta += generalisedCoordinates[rotOffset + 3*i + 1];
	}

	meanX = totalX/numberOfParticles;
	meanY = totalY/numberOfParticles;
	meanZ = totalZ/numberOfParticles;

>>>>>>> refs/remotes/origin/Michael
	meanAlpha = totalAlpha/numberOfParticles;
	meanBeta = totalBeta/numberOfParticles;

	// Calculate torques on each particle according to their alignment with average angle
	for (int i=0; i<numberOfParticles; i++)
	{
<<<<<<< HEAD
		// Calculate the distance between the particle and the average position
		dist = sqrt(pow(meanX - generalisedCoordinates[3*i + 0],2) + pow(meanY - generalisedCoordinates[3*i + 1],2) + pow(meanZ - generalisedCoordinates[3*i + 2],2));
		distMul = 1/dist; // Distance multiplier
		
		// Calculate torques in alpha and beta directions such that maximum torque is when the particle's axes are maximally separated from the average angle.
		additionalForces[rotOffset + 3*i + 0] += difAlpha = distMul * forceConst * sin(meanAlpha - generalisedCoordinates[rotOffset + 3*i + 0]);
		additionalForces[rotOffset + 3*i + 1] += difBeta = distMul * forceConst * sin(meanBeta - generalisedCoordinates[rotOffset + 3*i + 1]);
		/* Applies inverse distance multiplier as previously mentioned, which normalises the force so that it should be producing the average ~pi/2 angular change per timestep when the
			particle is in the order of a radius of the average position*/
		
		
		// Print stuff for debugging
		if (i<3) // Avoid some spam
		{
			printf("\ts:%e\tc:%e\tTheta:%e\tPhi:%e\n", dist, forceConst, generalisedCoordinates[rotOffset + 3*i + 0], generalisedCoordinates[rotOffset + 3*i + 1]);
			printf("Theta:\t%e\t%e\t%e\t%e\t%e\n", fmod(generalisedCoordinates[rotOffset + 3*i + 0], 2*gPi), distMul, forceConst, sin(meanAlpha - generalisedCoordinates[rotOffset + 3*i + 0]), additionalForces[rotOffset + 3*i + 0]);
			printf("Phi:\t%e\t%e\t%e\t%e\t%e\n", fmod(generalisedCoordinates[rotOffset + 3*i + 1], 2*gPi), distMul, forceConst, sin(meanBeta - generalisedCoordinates[rotOffset + 3*i + 1]), additionalForces[rotOffset + 3*i + 1]);
		}
=======
		distMul = 1/sqrt(pow(meanX - generalisedCoordinates[3*i + 0],2) + pow(meanY - generalisedCoordinates[3*i + 1],2) + pow(meanZ - generalisedCoordinates[3*i + 2],2));

		additionalForces[rotOffset + 3*i + 0] += difAlpha = distMul * forceConst * sin(meanAlpha - generalisedCoordinates[rotOffset + 3*i + 0]);
		additionalForces[rotOffset + 3*i + 1] += difBeta = distMul * forceConst * sin(meanBeta - generalisedCoordinates[rotOffset + 3*i + 1]);

		printf("%e\t%e\t%e\n", distMul, difAlpha, difBeta);
>>>>>>> refs/remotes/origin/Michael
	}

}

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField)
{

	for (int i=0; i<numberOfCells/6; i++)
	{
		additionalForces[3*i + 0] +=  drivingField.mag * abs(cos((drivingField.alpha - generalisedCoordinates[numberOfCells/2 + 3*i + 0])/2)) * abs(cos((drivingField.beta - generalisedCoordinates[numberOfCells/2 + 3*i + 1])/2));
		//additionalForces[3*i + 1] +=  drivingField.mag * cos(drivingField.beta - generalisedCoordinates[numberOfCells/2 + 3*i + 1]);
	}

}
