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
