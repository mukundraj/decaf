#ifndef PATHLINE_H
#define PATHLINE_H

#include <vector>
#include <Eigen/Dense>
#include <diy/master.hpp>
#include "ptrace.h"


class block;
class mpas_io;

class pathline
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

	void get_validated_cell_id(const block &mpas1, Eigen::Array3d &xSubStep, int &iCell, int &nCellVertices);

	void velocity_time_interpolation(const int timeInterpOrder, 
									const double *timeCoeff, 
									const int iCell, 
									const int iLevel, 
									const int nCellVertices, 
									const double zSubStep,
									block *b,
									mpas_io &mpas1, // includes verticesOnCell, boundaryVertex, zMid, zTop, xVertex, yVertex, zVertex
									const Eigen::Array3d &xSubStep,
									Eigen::Array3d &particleVelocity, 
									double &particleVelocityVert);

	void particle_vertical_treatment(const int nCellVertices, const int *verticesOnCell, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, block *b, mpas_io &mpas1,double uvCell[][3], const int iLevel, const double zLoc, const double *zMid, const double *zTop, const double *zVertVelocityTop, double &verticalVelocityInterp);

	Eigen::Array3d particle_horizontal_interpolation(mpas_io &mpas1, const int nCellVertices, const double vertCoords[][3], const Eigen::Array3d &pointInterp, const double uvCell[][3], double *areaB);

	double interp_vert_velocity_to_zlevel(block *b, mpas_io &mpas1, const int iLevel, const double zSubStep, const double *zTop, const double *vertVelocityTop);

	void interp_nodal_vectors(block *b, mpas_io &mpas1, const int nCellVertices, const int *verticesOnCell, const int iLevel, int nVertLevels, const double phiInterp, const double *phiVals, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, double uvCell[][3]);

	void get_bounding_indices(int &iHigh, int &iLow, const double phiInterp, const double *phiVals, const int iLevel, const int nVertLevels);
	void get_bounding_indices_brute_force(const int nVertLevels, const double phiInterp, const double *phiVals, int &iHigh, int &iLow);
	void wachspress_coordinates(const mpas_io &mpas1, const int nVertices, const double vertCoords[][3], const Eigen::Array3d &pointInterp, double *areaB, double *wachspress, int flag=0);

	double wachspress_interpolate(const double *lambda, const double phi[][3], const int idx, const int nCellVertices  );

	// computes the vertex velocites for all vertices that lie on cells in local_gcIds
	void cell_to_vertex_interpolation(block *b, std::vector<int> &local_gcIds);

	bool compute_streamlines(block *b, mpas_io &mpas1, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int prediction, size_t &nsteps, std::vector<EndPt> &particles_hold, int skip_rate);
	bool in_global_domain(const Pt& p, double cx, double cy, double cz);
	bool in_local_domain (const block *b, const Pt& p, int &iCell, int round);
};













#endif // PATHLINE_H


