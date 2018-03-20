
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "initial_finalisation.h"

extern int gDebug;
extern int gSerial;
extern int gNumOfthreads;
extern double gPi;


int cmd_line_read_in(int argc, char *argv[], environmentVariables *conditions)
{
	if (argc > 1)
	{
		for (int i=0; i<argc; i++)
		{
			if(strstr(argv[i],"-debug") != NULL) gDebug = 1;
			else if(strstr(argv[i],"-serial") != NULL)	gSerial = 1;

			else if (strstr(argv[i],"-num") != NULL || strstr(argv[i],"-n") != NULL)

			{
				if (sscanf(argv[i+1],"%d", &conditions->numberOfParticles) != 1)
				{
					printf("Invalid number of particles\n");
					return -1;
				}
			}
      		else if (strstr(argv[i],"-numthreads") != NULL)
			{
				if (sscanf(argv[i+1],"%d", &gNumOfthreads) != 1)
				{
					printf("Invalid number of threads\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-threads") != NULL)
			{
				if (sscanf(argv[i+1],"%d", &gNumOfthreads) != 1)
				{
					printf("Invalid number of threads\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-cube") != NULL)
			{
				double temp_num;
				if (strstr(argv[i+1],"R") != NULL)
				{
					if (sscanf(argv[i+1],"%lfR", &temp_num) != 2)
					{
						conditions->xMax = temp_num * conditions->radius;
						conditions->yMax = temp_num * conditions->radius;
						conditions->zMax = temp_num * conditions->radius;
					}
					else
					{
						printf("Invalid maximum dimension values\n");
						return -1;
					}
				}
				else if (sscanf(argv[i+1],"%lf", &temp_num) == 1)
				{
					conditions->xMax = temp_num;
					conditions->yMax = temp_num;
					conditions->zMax = temp_num;
				}
				else
				{
					printf("Invalid maximum dimension values\n");
					return -1;
				}

				printf("Set dimensions to %1.1em\n", conditions->xMax);

			}
			else if (strstr(argv[i],"-x") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->xMax) != 1)
				{
					printf("Invalid maximum x-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-y") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->yMax) != 1)
				{
					printf("Invalid maximum y-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-z") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->zMax) != 1)
				{
					printf("Invalid maximum z-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-t") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->endTime) != 1)
				{
					printf("Invalid duration\n");
					return -1;
				}
				else printf("Simulation duration set to %1.1e seconds\n", conditions->endTime);
			}
			else if (strstr(argv[i],"-dt") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->deltaTime) != 1)
				{
					printf("Invalid timestep\n");
					return -1;
				}
				else printf("Simulation timestep set to %1.1e seconds\n", conditions->deltaTime);
			}
		}
		if (gDebug == 1 && gSerial == 1) printf("Debug & serial modes active\n");
		else if (gDebug == 1 && gSerial == 0) printf("Debug mode active\n");
		else if (gDebug == 0 && gSerial == 1) printf("Serial mode active\n");

	}
	return 1;
}


void boilerplate_variables(environmentVariables *conditions)
{
    conditions->temperature = 298; // K
	conditions->viscosity = 8.9E-4; //N m^-2 s
	conditions->radius = 50E-9; // m
	conditions->currentTime = 0; // Seconds
	conditions->deltaTime = 1E-7; // Seconds
	conditions->endTime = 1E-2; // Seconds
	conditions->mass = (4/3) * gPi * pow(conditions->radius,3)*19320; // kg - density of gold
	conditions->xMax = 1E-7;
	conditions->yMax = 1E-7;
	conditions->zMax = 1E-7;
	conditions->numberOfParticles = 0;


	gNumOfthreads =omp_get_max_threads();
}
