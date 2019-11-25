#include "advect.h"
#include "mpas_io.h"
#include <cmath>
#include "misc.h"
#include "geometry_utils.h"
#include <algorithm>
#include <cstdlib>
#include "pathline.h"
#include <iostream>

// rk4 functions 



// mpas projection and associated function used by both streamline and pathline

void particle_horizontal_movement(const mpas_io& mpas1, Eigen::Array3d &pParticle, const Eigen::Array3d &dpParticle){


	double eps = 1e-10;
	double lenPath = sqrt((dpParticle*dpParticle).sum());

	double radiusShell = sqrt((pParticle*pParticle).sum());

	Eigen::Array3d pParticleTemp = pParticle + dpParticle;
	pParticleTemp = (radiusShell / sqrt((pParticleTemp * pParticleTemp).sum())) * pParticleTemp;

	double arcLen = get_arc_length(pParticle(0), pParticle(1), pParticle(2), 
		pParticleTemp(0), pParticleTemp(1), pParticleTemp(2));

	double alpha = 0;
	if (arcLen > eps){
		alpha = lenPath / arcLen;
		// dprint("alpha %f, arcLen %f", alpha, arcLen);
	}
	else return;


	// Eigen::Array3d pParticleInterp;
	spherical_linear_interp(pParticle, pParticle, pParticleTemp, alpha);

	// pParticle = pParticleInterp;


}



void spherical_linear_interp(Eigen::Array3d &pInterp, const Eigen::Array3d &p0, const Eigen::Array3d &p1, double alpha){

	double eps = 1e-14;

	double p0mag = sqrt((p0*p0).sum());
	double p1mag = sqrt((p1*p1).sum());

	Eigen::Array3d p0scaled = p0/p0mag;
	Eigen::Array3d p1scaled = p1/p1mag;
	std::cout.precision(15);
	// dprint("p0 (%.15f %.15f %.15f) p1 (%.15f %.15f %.15f)", p0scaled(0), p0scaled(1), p0scaled(2), p1scaled(0), p1scaled(1), p1scaled(2));

	double dotProd = std::min(1.0, std::max(-1.0, (p0scaled*p1scaled).sum()));

	double omega = acos(dotProd);

	
	if(std::abs(omega) < eps ){
		pInterp = p0;
		// dprint("omega %.15f, dotprod %.15f, eps %.15f", std::abs(omega), dotProd, eps);
		return;
	}

	pInterp = sin( (1-alpha) * omega) / sin(omega) * p0 + sin(alpha * omega) / sin(omega) * p1;
	// dprint("pInterp (%f %f %f)", pInterp(0), pInterp(1), pInterp(2));

}


