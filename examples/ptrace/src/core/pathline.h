#ifndef PATHLINE_H
#define PATHLINE_H

#include <vector>
#include <Eigen/Dense>
#include <diy/master.hpp>
#include "ptrace.h"
#include "flowline.h"


class block;
class mpas_io;

class pathline : public flowline
{

public:


	double dtSim, dtParticle, nmax;

	double nSteps, dt; // values set in constructor
	double radius;



	std::vector<double*> uVertexVelocities, vVertexVelocities, wVertexVelocities; // points to starting addresses of available velocities during current time step (initialize len=timeInterpOrder)
	double *zMid_cur, *zTop_cur, *vertVelocityTop_cur; // pointers to current time array start location of respective quantities

	/* For debugging */
	size_t zMid_cur_size, zTop_cur_size;

	std::vector<double> temp; // holder for zero values during the first advection round

	int timeInterpOrder = 2;

	pathline(mpas_io &mpas1_in, double dtSim_in, double dtParticle_in);

	void update_velocity_vectors(mpas_io &mpas1, int framenum);

	void compute_epoch(block *mpas1, int framenum);

	


	

	

	bool compute_flow(block *b, mpas_io &mpas1, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int prediction, size_t &nsteps, std::vector<EndPt> &particles_hold, int skip_rate);

};













#endif // PATHLINE_H


