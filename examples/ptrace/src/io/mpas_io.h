#ifndef MPAS_IO
#define MPAS_IO

#include <string>
#include <vector>
#include <map>
#include "nabo/nabo.h"
#include <diy/mpi.hpp>

class mpas_io
{

public:
	double radius; // radius of earth surface
	double rradius; // "real" radius: actual radius of a point
	double cradius; // radius of point coordinates on mesh
	double ratio;
	double cx, cy, cz;

	/* Field data */
	long long int nCells, nEdges, nVertices, nVertLevels, maxEdges, nVertLevelsP1;
    std::vector<double> latVertex, lonVertex, xVertex, yVertex, zVertex;
    std::vector<double> xyzCell, xCell, yCell, zCell;
    std::vector<int> indexToVertexID, indexToCellID;
    std::vector<int> verticesOnEdge, cellsOnVertex, verticesOnCell;

    std::vector<double> velocityX, velocityY, velocityZ;
    std::vector<std::vector<double>> velocityXv, velocityYv, velocityZv;
    std::vector<std::vector<double>> zTop, zMid, vertVelocityTop;
	std::vector<double> zTop_;
    std::vector<double> zTopVertex, zTopVertexNorm;
    std::map<int, int> vertexIndex, cellIndex;
    std::vector<int> nEdgesOnCell, maxLevelCell;
    std::vector<int> boundaryVertex;
    std::vector<int> cellsOnCell;

    // Global values
    size_t nVertices_global;

    /* LIGHT data */
	std::vector<double> xParticle, yParticle, zParticle;
	std::vector<double> xParticle_, yParticle_, zParticle_; // Mukund's tracer points
	size_t nParticles;
	size_t nTime;
	std::vector<double> zLevelParticle;


	mpas_io();

	void update_data(int data_id, int frame_no, std::vector<int> &data_int, std::vector<double> &data_dbl);
	void loadMeshFromNetCDF_CANGA(diy::mpi::communicator& world, const std::string& filename, long long int time_id);


	/* KDTree */
	Eigen::MatrixXd C_local;
	Nabo::NNSearchD*  nns_cells;
	void create_cells_kdtree();

};




#endif