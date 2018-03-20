#ifndef _INITIAL_FINALISATION_H
#define _INITIAL_FINALISATION_H

#include <stdlib.h>


typedef struct
{
	double temperature; // K
    double viscosity; //N m^-2 s
    double radius; // m
    double currentTime;
    double deltaTime; // Seconds
    double endTime; // Seconds
	double mass; // kg
    double xMax; //m
    double yMax; //m
    double zMax; //m
    int numberOfParticles;
}environmentVariables;

int cmd_line_read_in(int argc, char *argv[], environmentVariables *conditions);


#endif //_INITIAL_FINALISATION_H
