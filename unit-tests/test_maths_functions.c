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
