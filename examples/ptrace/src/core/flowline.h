#ifndef FLOWLINE_H
#define FLOWLINE_H

#include <Eigen/Dense>
#include "ptrace.h"
#include <diy/master.hpp>


class mpas_io;
class block;


class flowline{


public:

    void particle_vertical_treatment(const int nCellVertices, const int *verticesOnCell, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, block *b, double uvCell[][3], const int iLevel, const double zLoc, const double *zMid, const double *zTop, const double *zVertVelocityTop, double &verticalVelocityInterp);

	Eigen::Array3d particle_horizontal_interpolation(double radius, const int nCellVertices, const double vertCoords[][3], const Eigen::Array3d &pointInterp, const double uvCell[][3], double *areaB);

	double interp_vert_velocity_to_zlevel(block *b, const int iLevel, const double zSubStep, const double *zTop, const double *vertVelocityTop);

	void interp_nodal_vectors(block *b, const int nCellVertices, const int *verticesOnCell, const int iLevel, int nVertLevels, const double phiInterp, const double *phiVals, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, double uvCell[][3]);

	void get_bounding_indices(int &iHigh, int &iLow, const double phiInterp, const double *phiVals, const int iLevel, const int nVertLevels);
	void get_bounding_indices_brute_force(const int nVertLevels, const double phiInterp, const double *phiVals, int &iHigh, int &iLow);

    void get_validated_cell_id(const block &mpas1, Eigen::Array3d &xSubStep, int &iCell, int &nCellVertices);

	void velocity_time_interpolation(const int timeInterpOrder, 
									const double *timeCoeff, 
									const int iCell, 
									const int iLevel, 
									const int nCellVertices, 
									const double zSubStep,
									block *b,
									const Eigen::Array3d &xSubStep,
									Eigen::Array3d &particleVelocity, 
									double &particleVelocityVert);


	void cell_to_vertex_interpolation(block *b, bool pred_frame);

    void wachspress_coordinates( double radius, const int nVertices, const double vertCoords[][3], const Eigen::Array3d &pointInterp, double *areaB, double *wachspress, int flag=0);

	double wachspress_interpolate(const double *lambda, const double phi[][3], const int idx, const int nCellVertices  );

    // // computes the vertex velocites for all vertices that lie on cells in local_gcIds
	// void cell_to_vertex_interpolation(block *b, std::vector<int> &local_gcIds);

    bool in_global_domain(const Pt& p, double cx, double cy, double cz);
	

	virtual bool compute_flow(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int prediction, size_t &nsteps, std::vector<EndPt> &particles_hold, int skip_rate) = 0;

	virtual bool in_local_domain (const block *b, const Pt& p, int &iCell, int round, int pid) = 0;

};


#endif // FLOWLINE_H