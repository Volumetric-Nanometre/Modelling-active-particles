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
#include <stdio.h>
#include <stdlib.h>

#include "particles.h"

int main(int argc, char *argv[])
{
    int numberOfParticles=0;

    particleVariables* particles=NULL;
    //
    // Call function to read in particle data
    //
    if((numberOfParticles = particle_read_in(&particles))<=0)
    {
        return numberOfParticles;
    }


    //
    // Free memory
    //

    if(particles!=NULL)
    {
        free(particles);
        particles=NULL;
    }

    return 0;
}
