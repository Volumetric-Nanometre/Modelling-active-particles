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
#include "diffusionmatrix.h"
//

//

double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;


int main(int argc, char *argv[])
{
    int numberOfParticles = 0;

    particleVariables* particles = NULL;
    //
    // Call function to read in particle data
    //
    if( ( numberOfParticles = particle_read_in( &particles ) ) <= 0)
    {


        getchar();
        return numberOfParticles;
    }

    double *diffusionMatrix = NULL ;

    double temperature = 273; // Kelvin
    double viscosity = 1;
    double radius = 1;

    if( (diffusionMatrix = diffusion_matrix_creation( numberOfParticles, particles ,temperature, viscosity, radius)) == NULL )
    {
        free( particles );
        particles = NULL ;
        getchar();
        return errno;
    }

    printf("Diffusion matrix created\n" );
    getchar();

    //
    // Free memory
    //

    if( particles != NULL )
    {
        free( particles );
        particles = NULL;
    }
    if( diffusionMatrix != NULL )
    {
        free( diffusionMatrix );
        particles = NULL;
    }

    return 0;
}
