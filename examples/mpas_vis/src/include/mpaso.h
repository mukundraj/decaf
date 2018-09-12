#ifndef MPASO_H
#define MPASO_H

#include <string>
#include <map>
#include <vector>


class mpaso
{


	public:

		size_t nCells, nEdges, nVertices, nVertLevels;
		std::vector<double> latVertex, lonVertex, xVertex, yVertex, zVertex;
		std::vector<double> xyzCell; // xCell, yCell, zCell;
		std::vector<int> indexToVertexID, indexToCellID;
		std::vector<int> verticesOnEdge, cellsOnVertex;
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
};

#endif
