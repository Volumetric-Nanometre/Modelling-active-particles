// main file
#include <stdio.h>
#include <stdlib.h>

#include "particles.h"

int main()
{
    int numberOfParticles=0;

    particleVariables* particles=NULL;

    if((numberOfParticles = particle_read_in(&particles))<=0)
    {
        return numberOfParticles;
    }

    if(particles!=NULL)
    {
        free(particles);
        particles=NULL;
    }

    return 0;
}
