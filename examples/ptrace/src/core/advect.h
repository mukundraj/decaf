#ifndef ADVECT_H
#define ADVECT_H

#include <Eigen/Dense>

class mpas_io;



void particle_horizontal_movement(const mpas_io& mpas1, Eigen::Array3d &pParticle, const Eigen::Array3d &dpParticle);


void spherical_linear_interp(Eigen::Array3d &pInterp, const Eigen::Array3d &p0, const Eigen::Array3d &p1, double alpha);


#endif