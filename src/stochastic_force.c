
/*********************
* Date of creation 17/10/2017
* Author:
* Contact:
* Other Authors: N/A
**************************************
* History
**************************************/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "stochastic_force.h"


double *stochastic_displacement_creation(int numberOfParticles)
{
    double *stochastic_force = NULL ;
    //
    // Allocate and check memory
    //
    stochastic_force = calloc( 6 * numberOfParticles, sizeof ( *stochastic_force ));

    if( stochastic_force == NULL )
    {
        printf("-Error %d : %s\n", errno, strerror( errno ) );
        return NULL;
    }

    return stochastic_force;
}
