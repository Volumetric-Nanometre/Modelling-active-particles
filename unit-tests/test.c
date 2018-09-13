#include <stdio.h>
//
// files to be tested
//
#include "test_maths_functions.c"
#include <assert.h>


int main(int argv, char *argc[])
{
	printf("Unit tests for Active Particle Simulation.\n");
	printf("Unit tests will either mark as PASS or crash the test.\n");
	printf("------------------------------------------------------\n");
	printf("Mathematical Functions\n\n");
	test_kronecker_delta();
	test_levi_civita_density();
	test_linear_index_from_coordinates();
	test_coordinates_from_linear_index();
	return 0;
}
