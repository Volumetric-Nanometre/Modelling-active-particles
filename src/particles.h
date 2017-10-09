/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
*
* 09/10/2017 -Created particl_read_in prototype. Also created the struct
*             particleVariables to store the varibles in use. 
**************************************/


#ifndef _PARTICLES_H
#define _PARTICLES_H
#include <stdio.h>

typedef struct
{
    double x;
    double y;
    double z;
    double alpha;
    double beta;
    double gamma;

    double dx;
    double dy;
    double dz;
    double dalpha;
    double dbeta;
    double dgamma;
}particleVariables;


int  particle_read_in(particleVariables **particles);

#endif // _PARTICLES_H
