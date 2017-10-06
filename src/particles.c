#include "particles.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

int  particle_read_in(particleVariables **particles)
{
    int numberOfParticles=0;

    particleVariables* initParticles=NULL;

    //
    // Open and check input file
    //
    FILE* particleInput=fopen("particleInput.txt","r");

    if(particleInput == NULL)
    {
        printf("- Error %d : %s\n", errno,strerror(errno));
        numberOfParticles=-errno;
    }
    else
    {

        //
        // Get number of particles. If Illegal value, close file and return -1
        //

        if(fscanf(particleInput,"%d\n",&numberOfParticles)!=1 || numberOfParticles<=0 )
        {
            printf("- Error : Number of particles is invalid \n");
            numberOfParticles=0;
            fclose(particleInput);
        }
        else
        {
            //
            // Allocate memory for input then read until
            // # of particles reached or the last whole line is read.
            //

            if((initParticles=calloc(numberOfParticles,sizeof(*initParticles))) == NULL)
            {
                printf("- Error %d : %s\n", errno,strerror(errno));
                fclose (particleInput);
                numberOfParticles=-errno;
            }

            int i=0;

            while(i<numberOfParticles && fscanf(particleInput,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",&initParticles[i].x,&initParticles[i].y,&initParticles[i].z,&initParticles[i].dx,&initParticles[i].dy,&initParticles[i].dz,
                                                                &initParticles[i].alpha,&initParticles[i].beta,&initParticles[i].gamma,&initParticles[i].dalpha,&initParticles[i].dbeta,&initParticles[i].dgamma)==12)
            {
                i++;
            }

        }
    }
    fclose(particleInput);
    *particles=initParticles;
    return numberOfParticles;
}
