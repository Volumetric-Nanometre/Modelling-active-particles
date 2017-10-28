
#include "diffusionmatrix.h"
#include "maths_functions.h"
#include "particles.h"

extern double gBoltzmannConst; // m^2 kg s^-2 K^-1
extern int gDebug;

void initialVelocities(int numberOfParticles, particleVariables *particles, environmentVariables *conditions, long int tSeed)
{

	double totalEnergy = gBoltzmannConst * conditions->temperature;

	if (gDebug) printf("Total energy of each particle is %.2eJ\n",totalEnergy);

	// kT = (1/2)mv^2 + (1/2)Iw^2
	// => v^2 = kT/m  && w^2 = kT/I		- if split equally between rotational and translational kinetic energies



	// randomly generate components given these constraints:
	// v^2 = v_x^2 + v_y^2 + v_z^2
	// w^2 = w_x^2 + w_y^2 + w_z^2

	for (int i = 0; i < numberOfParticles; i++)
	{

		double frac = ran1(&tSeed); // fraction of kinetic energy which is translational
		double vSquared = totalEnergy*frac/(0.5 * conditions->mass);
		double wSquared = totalEnergy*(1-frac)/( 0.5 * 0.4 * conditions->mass * pow(conditions->radius,2) ); // Moment of inertia of a solid sphere is (2/5)mr^2

		particles[i].dx = randSign(tSeed) * sqrt( vSquared ) * ran1(&tSeed); // randomly generate v_x^2 as a fraction of v^2
		particles[i].dy = randSign(tSeed) * sqrt( vSquared - pow(particles[i].dx,2) ) * ran1(&tSeed); // randomly generate v_y^2 as a fraction of v^2 - v_x^2
		particles[i].dz = randSign(tSeed) * sqrt( vSquared - pow(particles[i].dx,2) - pow(particles[i].dy,2) ); // calculate v_z^2 so that v_x^2 + v_y^2 + v_z^2 = v^2

		particles[i].dalpha = randSign(tSeed) * sqrt( wSquared ) * ran1(&tSeed); // randomly generate w_x^2 as a fraction of r^2
		particles[i].dbeta = randSign(tSeed) * sqrt( wSquared - pow(particles[i].dalpha,2) ) * ran1(&tSeed); // randomly generate w_y^2 as a fraction of w^2 - w_x^2
		particles[i].dgamma = randSign(tSeed) * sqrt( wSquared - pow(particles[i].dalpha,2) - pow(particles[i].dbeta,2) ); // calculate w_z^2 so that w_x^2 + w_y^2 + w_z^2 = w^2

		if (gDebug && i<3)
		{
			printf("Particle %d: v=(%.2e, %.2e, %.2e) w=(%.2e, %.2e, %.2e)\n", i, particles[i].dx, particles[i].dy, particles[i].dz, particles[i].dalpha, particles[i].dbeta, particles[i].dgamma);
			double TKE = 0.5 * conditions->mass * (pow(particles[i].dx,2) + pow(particles[i].dy,2) + pow(particles[i].dz,2));
			double RKE = 0.5 * 0.4 * conditions->mass * pow(conditions->radius,2) * (pow(particles[i].dalpha,2) + pow(particles[i].dbeta,2) + pow(particles[i].dgamma,2));
			printf("Trans KE: %.2eJ, rot KE: %.2eJ, total KE: %.2eJ\n",TKE,RKE,TKE+RKE);
		}
	}


}
