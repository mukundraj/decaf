#ifndef MPASO_H
#define MPASO_H

#include <string>
#include <map>
#include <vector>
#include <set>

#include "nabo/nabo.h"


class mpaso
{


	public:

		size_t nCells, nEdges, nVertices, nVertLevels;
		std::map<int, int> cgid_to_bid;

		std::vector<double> latVertex, lonVertex, xVertex, yVertex, zVertex;
		std::vector<double> xyzCell; // xCell, yCell, zCell;
		std::vector<int> indexToVertexID, indexToCellID;
		std::vector<int> verticesOnEdge, cellsOnVertex, verticesOnCell;
		std::vector<double> velocityX, velocityY, velocityZ;
		std::vector<double> velocityXv, velocityYv, velocityZv;
		std::vector<double> zTop;
		std::vector<double> zTopVertex, zTopVertexNorm;
		std::map<int, int> vertexIndex, cellIndex;

		std::vector<double> xCells, yCells, zCells;


		std::vector<std::vector<int> > dd_adjmat; 
		mpaso();
		~mpaso();
		void generate_domain_decomposition_graph(std::string &filename, int nblocks);
		void loadMeshFromNetCDF_CANGA(const std::string& filename, size_t time_id=0);
		void load_mesh_from_decaf_data_bar(std::vector<double> &decaf_block);
		void compute_cellIndex(int cblock, int nblocks);
		void read_cell_g_neighbors(); // reads graph.info and populates cell_g_neighbors
		

		int get_bid_for_pt(double *coords);


  		Nabo::NNSearchD*  nns_cells;
		Nabo::NNSearchD* nns;

		//Nabo::NNSearchD* nns_cells_local;
		//Nabo::NNSearchD* nns_local;



		std::vector<int> indexToCellID_p;
		std::vector<double> xCell_p;
		std::vector<double> yCell_p;
		std::vector<double> zCell_p;

		std::vector<double> velocityX_p;
		std::vector<double> velocityY_p;
		std::vector<double> velocityZ_p;

		std::vector<int> indexToVertexID_p;
		std::vector<double> xVertex_p;
		std::vector<double> yVertex_p;
		std::vector<double> zVertex_p;
		std::vector<double> zTop_p;
		std::vector<int> cellsOnVertex_p; 

		std::vector<std::vector<int>> cell_g_neighbors;

		std::set<int> ghost_cell_ids;
		std::set<int> block_vert_ids; // global ids of vertices native to block


		int init, done;

		int gid;
		
		double radius;	
};

#endif
