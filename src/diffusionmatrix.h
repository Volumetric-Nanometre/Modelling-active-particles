#ifndef _DIFFUSIONMATRIX_H
#define _DIFFUSIONMATRIX_H

#include <stdio.h>
#include "particles.h"

double *diffusion_matrix_creation(int numberOfParticles, particleVariables *particles, double temperature, double viscosity, double radius);

void translational_tensor_creation(double *oseenMatrix, particleVariables *particles, double temperature, double viscosity, double radius, int i, int j);

#endif //_DIFFUSIONMATRIX_H
