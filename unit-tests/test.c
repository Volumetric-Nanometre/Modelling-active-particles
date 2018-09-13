#include <stdio.h>
//
// files to be tested
//
#include "test_maths_functions.c"
#include "test_forces.c"
#include <assert.h>

double gBoltzmannConst = 1-23; // m^2 kg s^-2 K^-1
double gPi = 3;
double gGrav = 10; // m s^-2
int gDebug = 0;

int main(int argv, char *argc[])
{
	printf("Unit tests for Active Particle Simulation.\n");
	printf("Unit tests will either mark as PASS or crash the test.\n");
	printf("------------------------------------------------------\n\n");
	printf("Mathematical Functions\n\n");
	test_kronecker_delta();
	test_levi_civita_density();
	test_linear_index_from_coordinates();
	test_coordinates_from_linear_index();

	printf("\nForce Functions\n\n");
	test_force_gravity();
	return 0;
}
