#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#include "../src/forces.c"
//
// Test the force of gravity is correct when g=10
//

void test_force_gravity()
{
	printf("%s","force_gravity..." );

	double testForces[6]={1,2,3,4,5,6};


	force_gravity(testForces, 6, 1);

	assert(testForces[0] == 1); //x
	assert(testForces[1] == 2); //y
	assert(testForces[2] == -7); //z - gravity should only appear here. current value - mass*grav = 3 - 1*10 = -7
	assert(testForces[3] == 4); //alpha
	assert(testForces[4] == 5); //beta
	assert(testForces[5] == 6); //gamma
	printf("[PASS]\n");
}

/*void test_force_van_der_waals()
{
	printf("%s","force_van_der_waals..." );

	double testForces[12]={1,2,3,4,5,6,7,8,9,10,11,12};
	double coordinates[12]={1,2,3,2,4,6,0,0,0,0,0,0};

	force_van_der_waals(testForces, coordinates, 12, 1);

	assert(testForces[0] == 1); //x
	assert(testForces[1] == 2); //y
	assert(testForces[2] == -7); //z - gravity should only appear here. current value - mass*grav = 3 - 1*10 = -7
	assert(testForces[3] == 4); //alpha
	assert(testForces[4] == 5); //beta
	assert(testForces[5] == 6); //gamma
	printf("[PASS]\n");
}*/
/*
static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double radius);

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, int numberOfCells);

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells, environmentVariables conditions);

static void viseck_alignment_torque(double *additionalForces, double *generalisedCoordinates, environmentVariables conditions, double rCutoff);

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField);

static void polar_driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double drivingForceMagnitude);
*/
