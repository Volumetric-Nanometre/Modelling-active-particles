/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
*
**************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "maths_functions.h"

//
// Performs basic kronecker delta
//

int kronecker_delta(int i, int j)
{
    if(i == j)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

//
// Performs the levi_civita_density permutation
//

int levi_civita_density(int i, int j)
{
    if(i == j)
    {
        return 0;
    }
    else if( i < j)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}
