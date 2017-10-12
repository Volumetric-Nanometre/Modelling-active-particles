/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* Change History
**************************************/
#include <stdio.h>
#include <stdlib.h>

#include "particles.h"
#include "diffusionmatrix.h"
//

//

double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;

int gDebug = 0;

int main(int argc, char *argv[])
{
    //
    // Check Debug mode
    //

    //if(argc >= 2)
    {
        gDebug = 1;
        printf("Warning debug mode entered. Press any key to continue...\n" );
        getchar();
    }


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

    printf("Data read in success\n" );

    double *diffusionMatrix = NULL ;

    double temperature = 273; // Kelvin
    double viscosity = 1;
    double radius = 1;

    //
    // Create diffusion matrix
    //

    if( (diffusionMatrix = diffusion_matrix_creation( numberOfParticles, particles ,temperature, viscosity, radius)) == NULL )
    {
        free( particles );
        particles = NULL ;
        getchar();
        return errno;
    }

    printf("Diffusion matrix created\n" );

    //---------------------------- DEBUG------------------------------//
    //
    // Prints the diffusionMatrix to a file for inspection
    //
    if( gDebug == 1)
    {
        FILE *matrixOutput = fopen( "matrix_output.txt","w");

        for(int i = 0; i < 6 * numberOfParticles; i++)
        {
            for(int j = 0; j < 6 * numberOfParticles; j++)
            {
                fprintf(matrixOutput, "%e\t", diffusionMatrix[i * 6 * numberOfParticles + j]);
            }
            fprintf(matrixOutput, "\n");

        }
    }
    //---------------------------END---------------------------------//



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
