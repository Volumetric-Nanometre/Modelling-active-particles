/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: Oliver Hinds
**************************************
* History
*
**************************************/

#ifndef _MATHS_FUNCTIONS_H
#define _MATHS_FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int kronecker_delta(int i, int k);

int levi_civita_density(int i, int j);

double randSign(time_t tSeed);

float ran1(long *idum);

float guassdev(long *idum);


#endif  // _MATHS_FUNCTIONS_H
