#ifndef _DIFFUSIONMATRIX_H
#define _DIFFUSIONMATRIX_H

#include <stdio.h>

double *diffusion_matrix_creation(int numberOfParticles, double temperature, double viscosity, double radius);

void oseen_tensor_creation(double *oseenMatrix, double temperature, double viscosity, double radius, int i, int j);

#endif //_DIFFUSIONMATRIX_H
