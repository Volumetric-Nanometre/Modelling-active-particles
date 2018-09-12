#include <stdio.h>
//
// files to be tested
//
#include "test_maths_functions.c"
#include <assert.h>


int main(int argv, char *argc[])
{
	test_kronecker_delta();
	test_levi_civita_density();
	return 0;
}
