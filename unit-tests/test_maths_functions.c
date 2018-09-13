#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#include "../src/maths_functions.c"
//
// Test the Kronecker delta outputs 1 when a==b and 0 when a!=b
//

void test_kronecker_delta()
{
	printf("%s","kronecker_delta..." );

	assert(kronecker_delta(1, 1)==1);
	assert(kronecker_delta(1, 0)==0);
	assert(kronecker_delta(0, 1)==0);
	printf("[PASS]\n");
}

//
// Test the levi_civita_density(int i, int j) delta outputs 1 when a==b and 0 when a!=b
//

void test_levi_civita_density()
{
	printf("%s","levi_civita_density..." );

	assert(levi_civita_density(0, 1)==1);
	assert(levi_civita_density( 1,  0)==-1);
	assert(levi_civita_density( 2,  1)==-1);
	assert(levi_civita_density( 2,  2)==0);
	printf("[PASS]\n");
}

//
// Test the linear_index_from_coordinates() function. x = 3 y = 2 z = 1. Therefore
// 25+10+3 = 38
//

void test_linear_index_from_coordinates()
{
	printf("%s","linear_index_from_coordinates..." );

	assert(linear_index_from_coordinates(5, 5, 3,  2,  1) == 38);

	printf("[PASS]\n");
}

//
// Test the test_coordinates_from_linear_index() function.Index of 38, therefore x = 3 y = 2 z = 1.
//

void test_coordinates_from_linear_index()
{
	int x=0, y=0, z=0;
	printf("%s","coordinates_from_linear_index..." );

	coordinates_from_linear_index(38, 5, 5, &x, &y, &z);

	assert(x == 3 && y == 2 && z == 1);

	printf("[PASS]\n");

}
