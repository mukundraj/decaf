#include "mpaso.h"
#include <netcdf.h>
#include "def.h"
#include "interpolators.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "misc.h"

#include <Eigen/Dense>

mpaso::mpaso(){
	radius = 6371220.;
	init = 0;
	done = 0;
}

mpaso::~mpaso(){
}

// int mpaso::get_bid_for_pt(double *coords){
	
//  	int bid;	
	
// 	// find nearest cell
// 	Eigen::VectorXi nearest_cell_idx(1);
// 	Eigen::VectorXd dists2(1);
// 	Eigen::VectorXd q(3);

// 	q << coords[0], coords[1] , coords[2];		
// 	nns_cells->knn(q, nearest_cell_idx,dists2, 1);
// 	// dprint("nearest_cell_idx cgid %d %d %d %d", nearest_cell_idx[0], indexToCellID[nearest_cell_idx[0]], cgid_to_bid[indexToCellID[nearest_cell_idx[0]]], cgid_to_bid[4463]);
// 	// find the bid of the cell
// 	bid = cgid_to_bid[indexToCellID[nearest_cell_idx[0]]];
// 	return bid;		
// }

// TODO: check if following function is needed. Why not use indexToCellID like for vertices?
// mapping from global cell id to local cell id
void mpaso::compute_cellIndex(int cblock, int nblocks){

	// now populate the cellID_to_bid vector
	std::string ip_file = "graph.info.part."+itos(nblocks);
	std::ifstream file;
	file.clear();
	
	int g_idx=0, cblock_idx=0, temp;
	file.open(ip_file);

	std::string line;

	while (std::getline(file,line)) { 
		std::istringstream iss(line);
		iss >>  temp;
		// fprintf(stderr, " temp %d %d\n", temp, cblock);
		if (temp==cblock){
			cellIndex[g_idx+1] = cblock_idx;
			cblock_idx++;
		}


		g_idx++; 
	}

	file.close();


}

void mpaso::read_cell_g_neighbors(){

	std::string ip_file = "graph.info";

	std::ifstream file(ip_file.c_str());
	// file.clear();
	std::string line;
	int temp;
	std::getline(file,line); // reading the header line
	std::istringstream iss(line);
	iss >>  temp;
	cell_g_neighbors.resize(temp+1); 


	int cur_cell_1idx = 1;
	while (std::getline(file,line)) { 
		std::istringstream iss2(line);
		while(iss2 >>  temp){

			cell_g_neighbors[cur_cell_1idx].push_back(temp);
		}
		cur_cell_1idx++;
		
		// fprintf(stderr, " cur_cell_1idx %d\n", cur_cell_1idx);
		
	}

	// for (int i=0;i<cell_g_neighbors.size();i++){
	// 	fprintf(stderr, "%d %ld\n", i, cell_g_neighbors[i].size());
	// }
	// fprintf(stderr, "%d %d %d %d %d %d\n", cell_g_neighbors[7234][0], cell_g_neighbors[7234][1], 
	// 					cell_g_neighbors[7234][2], cell_g_neighbors[7234][3], cell_g_neighbors[7234][4], cell_g_neighbors[7234][5]);


}

void mpaso::generate_domain_decomposition_graph(std::string &filename, int nblocks){

	// fprintf(stderr, "%s \n", filename.c_str());
	// std::ifstream file(filename.c_str());
	
	// std::string line;
 //        // Read one line at a time into the variable line:
 //        while(std::getline(file, line))
	// {
	// 	std::vector<int>   lineData;
	// 	std::stringstream  lineStream(line);

	// 	int value;
	// 	// Read an integer at a time from the line
	// 	while(lineStream >> value)
	// 	{
	// 		// Add the integers from a line to a 1D array (vector)
	// 		lineData.push_back(value);

	// 	}
	// 	this->dd_adjmat.push_back(lineData);
	// }
	// file.close();
	// fprintf(stderr, "dd graph %d %d %d %d", dd_adjmat[0][0], dd_adjmat[0][1], dd_adjmat[1][0], dd_adjmat[1][1]);

	// file.clear();
	// now populate the cellID_to_bid vector
	std::string ip_file = "graph.info.part."+itos(nblocks);
	int ctr=0, temp;
	// file.open(ip_file);
	// while ( !file.eof ()  ) {   
	// 	ctr++; 
	// 	file >> temp;
	// }
	// file.close();
	// cgid_to_bid.resize(ctr);
	std::ifstream file;
	file.open(ip_file);
	std::string line;
	ctr=1;
	while (std::getline(file,line)) { 
		std::istringstream iss(line);
		iss >>  temp;
		cgid_to_bid[ctr] = temp;
		if (ctr==4463)
			 dprint("ctr cgid_to_bid[i] %d %d", ctr, cgid_to_bid[ctr]);
		ctr++;
	}
	file.close();

}
void mpaso::load_mesh_from_decaf_data_bar(std::vector<double> &data_bar){

	//std::cout<<"inside new function "<<"\n";
	//printf("fir cell vals %f\n", data_bar[int(12+data_bar[4]+data_bar[5]+data_bar[6]+data_bar[7])]);

	int n_procs = int(data_bar[17]);
	int frame_no = int(data_bar[13]); 
	std::vector<int> bar_start_inds(n_procs);
	nVertices = 0;
	nCells = 0;
	nVertLevels= 0;

	//std::vector<double> xCells, yCells, zCells; 


	bar_start_inds[0] = 0;
	// get bar start positions
	for (int i=0; i< n_procs; i++){
		//nVertices += int(data_bar[20+i*]);
		if (i>0){
			bar_start_inds[i] = bar_start_inds[i-1] + int(data_bar[bar_start_inds[i-1]+18]);
		}
		//std::cout<<"mpas frame no: "<<data_bar[bar_start_inds[i]+13]<<"\n";
	}

	// get the total lengths
	for (int i=0; i<n_procs; i++){
		nVertices += int(data_bar[int(bar_start_inds[i]+14)]);
		nCells += int(data_bar[int(bar_start_inds[i]+15)]);
		nVertLevels = int(data_bar[int(bar_start_inds[i]+16)]);
	}

	//std::cout<<"nVertices "<<nVertices<<" "<<nCells<<" "<<nVertLevels<< "\n";   
	// resize the arrays

	indexToCellID.resize(nCells);	 
	velocityX.resize(nVertLevels*nCells);
	velocityY.resize(nVertLevels*nCells);
	velocityZ.resize(nVertLevels*nCells);
	zTop.resize(nVertLevels*nCells);

	//	if (cellsOnVertex.size()==0){
	//cellsOnVertex.resize(nVertices*3);
	//}
	int read_order[13] = {4, 5, 6, 11, 0, 1, 2, 3, 8, 9,10,7, 12};
	int field_offsets[13];
	int field_sizes[13];


	// two loops: one for dynamic and another for static, read in order of write
	// loop through procs and fill arrays
	for (int i=0;i<n_procs;i++){
		int bsp = bar_start_inds[i]; //bar start position for current proc
		//std::cout<<"mpaso initial bar values "<<data_bar[bsp+0]<<" "<<data_bar[bsp+4]<<"\n";
		// collect field global offsets and sizes in bar
		field_offsets[read_order[0]] = bsp+20;
		field_sizes[read_order[0]] = data_bar[bsp+read_order[0]];
		for (int j=1; j<13; j++){
			field_offsets[read_order[j]] = field_offsets[read_order[j-1]]+field_sizes[read_order[j-1]];
			field_sizes[read_order[j]] = data_bar[bsp+read_order[j]];
		}


		// if frame==1, parse through index to cellID 
		// if frame==1 parse through vertex order for cellsOnVertex   
		if (frame_no==1){
		//	std::cout<<"populating frame 1  values "<<frame_no<<"\n";
			cellsOnVertex.resize(nVertices*3);
			// verticesOnCell.resize(nCells*6);
			xCells.resize(nCells);
			yCells.resize(nCells);
			zCells.resize(nCells);
			indexToVertexID.resize(nVertices);


			xVertex.resize(nVertices);
			yVertex.resize(nVertices);
			zVertex.resize(nVertices);

			//for (int j=0; j<13; j++){
			//std::cout<<j<<" "<<read_order[j]<<" "<<field_offsets[read_order[j]]<<" "<<field_sizes[read_order[j]]<<"\n";
			// }
			int k = 0; // xCell index offset var
			for (int j=field_offsets[0]; j<field_offsets[0]+field_sizes[0]; j++){
				xCells[data_bar[j]-1] = data_bar[field_offsets[1]+k]; // -1 on LHS to deal with fortran to C offset
				yCells[data_bar[j]-1] = data_bar[field_offsets[2]+k];
				zCells[data_bar[j]-1] = data_bar[field_offsets[3]+k];
				indexToVertexID[data_bar[j]-1] = data_bar[j]; 
				//std::cout<<" "<<int(data_bar[j]);
				k++;
			}



			k = 0;
			for (int j=field_offsets[7];j<field_offsets[7]+field_sizes[7];j++){
				xVertex[data_bar[j]-1] = data_bar[field_offsets[8]+k];
				yVertex[data_bar[j]-1] = data_bar[field_offsets[9]+k];
				zVertex[data_bar[j]-1] = data_bar[field_offsets[10]+k];
				k++;
			}

			//	for (int j=0;j<100; j++){
			//		std::cout<<xVertex[j]<<"\n ";
			//	}
			int l=0;
			for (int j=field_offsets[7];j<field_offsets[7]+field_sizes[7];j++){
				int m = data_bar[j] - 1;
				for (int k=0; k<3; k++){
					cellsOnVertex[m*3+k] = data_bar[field_offsets[12]+3*l+k];
				}
				l++;	
			}

			

		}

		/*	if (frame_no==2){
			for (int j=0; j<13; j++){
			std::cout<<j<<" "<<read_order[j]<<" "<<field_offsets[read_order[j]]<<" "<<field_sizes[read_order[j]]<<"\n";
			}
			}
		 */

		//std::cout<<"populating non frame 1 related values\n";
		int l=0;
		for (int j=field_offsets[0]; j<field_offsets[0]+field_sizes[0];j++){
			int m = data_bar[j]-1;
			for (int k=0;k<nVertLevels;k++){
				velocityX[m*nVertLevels+k] = data_bar[field_offsets[4]+nVertLevels*l+k];
				velocityY[m*nVertLevels+k] = data_bar[field_offsets[5]+nVertLevels*l+k];	
				velocityZ[m*nVertLevels+k] = data_bar[field_offsets[6]+nVertLevels*l+k];
				zTop[m*nVertLevels+k] = data_bar[field_offsets[11]+nVertLevels*l+k];
				indexToCellID[m] = data_bar[j];
			}
			l++;
		}




	}

	//std::cout<<"zCell size "<<frame_no<<" "<<zCells.size()<<" ";
	//for (int j=0;j<110;j++){
	//	std::cout<<indexToCellID[j]<<" ";
	//}

	// loop through procs and fill vertex based arrays
	//for (int i=0; i< 230; i++){
	//std::cout<<i<<" "<<data_ar[bar_start_inds[0]+i]<<"\n";
	//}

	// now populating derived fields needed for flow visualization
	//vertexIndex, cellIndex,  zTopVertex, zTopVertexNorm, cellIndex, velocityXv:

	for (int i=0; i<nCells; i++) { // shoudn't i vary from 1 to nCells here? then cellIndex[indexToCellID[i-1]] = i-1
		cellIndex[indexToCellID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
	//	std::cout<<i<<" "<<indexToCellID[i]<<" "<<cellIndex[i]<<"\n";
	}
	for (int i=0; i<nVertices; i++) {
		vertexIndex[indexToVertexID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	}

	// derive velocity on cell verticies (on only one altitude level here)
	velocityXv.resize(nVertices*nVertLevels);
	velocityYv.resize(nVertices*nVertLevels);
	velocityZv.resize(nVertices*nVertLevels);
	zTopVertex.resize(nVertices*nVertLevels);
	zTopVertexNorm.resize(nVertices*nVertLevels);
	xyzCell.resize(nCells*3);
	for (int i=0; i<nCells; i++) xyzCell[i*3] = xCells[i];
	for (int i=0; i<nCells; i++) xyzCell[i*3+1] = yCells[i];
	for (int i=0; i<nCells; i++) xyzCell[i*3+2] = zCells[i];

	for (int i=0; i<nVertices; i++) {
		//    if (cellsOnVertex[i*3] == 0 || cellsOnVertex[i*3+1] == 0 || cellsOnVertex[i*3+2] == 0) continue; // on boundary
		int c0 = cellIndex[cellsOnVertex[i*3]], c1 = cellIndex[cellsOnVertex[i*3+1]], c2 = cellIndex[cellsOnVertex[i*3+2]];

		double X[3][3] = {
			{xyzCell[c0*3], xyzCell[c0*3+1], xyzCell[c0*3+2]},
			{xyzCell[c1*3], xyzCell[c1*3+1], xyzCell[c1*3+2]},
			{xyzCell[c2*3], xyzCell[c2*3+1], xyzCell[c2*3+2]},

		};

		double P[3] = {xVertex[i], yVertex[i], zVertex[i]};

		double lambda[3];
		barycentric_point2triangle(X[0], X[1], X[2], P, lambda);
		//[i * nVertLevels + curVertLevel];
		for (int curVertLevel=0;curVertLevel<nVertLevels; curVertLevel++){
		
			velocityXv[i*nVertLevels+curVertLevel] = lambda[0] * velocityX[c0*nVertLevels+curVertLevel] + lambda[1] * velocityX[c1*nVertLevels+curVertLevel] + lambda[2] * velocityX[c2*nVertLevels+curVertLevel];
			velocityYv[i*nVertLevels+curVertLevel] = lambda[0] * velocityY[c0*nVertLevels+curVertLevel] + lambda[1] * velocityY[c1*nVertLevels+curVertLevel] + lambda[2] * velocityY[c2*nVertLevels+curVertLevel];
			velocityZv[i*nVertLevels+curVertLevel] = lambda[0] * velocityZ[c0*nVertLevels+curVertLevel] + lambda[1] * velocityZ[c1*nVertLevels+curVertLevel] + lambda[2] * velocityZ[c2*nVertLevels+curVertLevel];

			zTopVertex[nVertLevels*i + curVertLevel] = lambda[0] * zTop[c0*nVertLevels+curVertLevel] + lambda[1] * zTop[c1*nVertLevels+curVertLevel] + lambda[2] * zTop[c2*nVertLevels+curVertLevel];


		}


	}

printf("nVertLevels %ld\n", nVertLevels);
	// then streamlines

	// then a manager function that computes the num of steps and calls pathline filter
}

void mpaso::loadMeshFromNetCDF_CANGA(const std::string& filename, size_t time_id){


	int ncid;
	int dimid_cells, dimid_edges, dimid_vertices, dimid_vertLevels;
	int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex,
	    varid_latCell, varid_lonCell, varid_xCell, varid_yCell, varid_zCell,
	    varid_verticesOnEdge, varid_cellsOnVertex,
	    varid_indexToVertexID, varid_indexToCellID,
	    varid_velocityX, varid_velocityY, varid_velocityZ;
	int varid_zTop;

	NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

	NC_SAFE_CALL( nc_inq_dimid(ncid, "nCells", &dimid_cells) );
	NC_SAFE_CALL( nc_inq_dimid(ncid, "nEdges", &dimid_edges) );
	NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_vertices) );
	NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertLevels", &dimid_vertLevels) );

	NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_cells, &nCells) );
	NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_edges, &nEdges) );
	NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertices, &nVertices) );
	NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertLevels, &nVertLevels) );

	NC_SAFE_CALL( nc_inq_varid(ncid, "indexToVertexID", &varid_indexToVertexID) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "indexToCellID", &varid_indexToCellID) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "latCell", &varid_latCell) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "lonCell", &varid_lonCell) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "xCell", &varid_xCell) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "yCell", &varid_yCell) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "zCell", &varid_zCell) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "latVertex", &varid_latVertex) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "lonVertex", &varid_lonVertex) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "verticesOnEdge", &varid_verticesOnEdge) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "cellsOnVertex", &varid_cellsOnVertex) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "velocityX", &varid_velocityX) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "velocityY", &varid_velocityY) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "velocityZ", &varid_velocityZ) );
	NC_SAFE_CALL( nc_inq_varid(ncid, "zTop", &varid_zTop));

	const size_t start_cells[1] = {0}, size_cells[1] = {nCells};

	indexToCellID.resize(nCells);
	NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToCellID, start_cells, size_cells, &indexToCellID[0]) );
	for (int i=0; i<nCells; i++) {
		cellIndex[indexToCellID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
	}

	std::vector<double> coord_cells;
	coord_cells.resize(nCells);
	xyzCell.resize(nCells*3);
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &coord_cells[0]) );
	for (int i=0; i<nCells; i++) xyzCell[i*3] = coord_cells[i];
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &coord_cells[0]) );
	for (int i=0; i<nCells; i++) xyzCell[i*3+1] = coord_cells[i];
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &coord_cells[0]) );
	for (int i=0; i<nCells; i++) xyzCell[i*3+2] = coord_cells[i];

	const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
	latVertex.resize(nVertices);
	lonVertex.resize(nVertices);
	xVertex.resize(nVertices);
	yVertex.resize(nVertices);
	zVertex.resize(nVertices);
	indexToVertexID.resize(nVertices);

	NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToVertexID, start_vertices, size_vertices, &indexToVertexID[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_latVertex, start_vertices, size_vertices, &latVertex[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_lonVertex, start_vertices, size_vertices, &lonVertex[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, &xVertex[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, &yVertex[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, &zVertex[0]) );

	for (int i=0; i<nVertices; i++) {
		vertexIndex[indexToVertexID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	}

	const size_t start_edges2[2] = {0, 0}, size_edges2[2] = {nEdges, 2};
	verticesOnEdge.resize(nEdges*2);

	NC_SAFE_CALL( nc_get_vara_int(ncid, varid_verticesOnEdge, start_edges2, size_edges2, &verticesOnEdge[0]) );

	// for (int i=0; i<nEdges; i++)
	//   fprintf(stderr, "%d, %d\n", verticesOnEdge[i*2], verticesOnEdge[i*2+1]);

	const size_t start_vertex_cell[2] = {0, 0}, size_vertex_cell[2] = {nVertices, 3};
	cellsOnVertex.resize(nVertices*3);

	NC_SAFE_CALL( nc_get_vara_int(ncid, varid_cellsOnVertex, start_vertex_cell, size_vertex_cell, &cellsOnVertex[0]) );

	// getting velocity at one time at cell centers
	const size_t start_time_cell_level[3] = {time_id, 0, 0}, size_time_cell_level[3] = {1, nCells, nVertLevels};
	velocityX.resize(nCells*nVertLevels);
	velocityY.resize(nCells*nVertLevels);
	velocityZ.resize(nCells*nVertLevels);



	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityX, start_time_cell_level, size_time_cell_level, &velocityX[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityY, start_time_cell_level, size_time_cell_level, &velocityY[0]) );
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityZ, start_time_cell_level, size_time_cell_level, &velocityZ[0]) );

	zTop.resize(nCells*nVertLevels);
	zTopVertex.resize(nVertices*nVertLevels);
	zTopVertexNorm.resize(nVertices*nVertLevels);

	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zTop, start_time_cell_level, size_time_cell_level, &zTop[0]) );

	// derive velocity on cell verticies (on only one altitude level here)
	velocityXv.resize(nVertices*nVertLevels);
	velocityYv.resize(nVertices*nVertLevels);
	velocityZv.resize(nVertices*nVertLevels);

	//  std::cout<<time+1<<" time here";
	//  exit(0);

	for (int i=0; i<nVertices; i++) {
		//    if (cellsOnVertex[i*3] == 0 || cellsOnVertex[i*3+1] == 0 || cellsOnVertex[i*3+2] == 0) continue; // on boundary
		int c0 = cellIndex[cellsOnVertex[i*3]], c1 = cellIndex[cellsOnVertex[i*3+1]], c2 = cellIndex[cellsOnVertex[i*3+2]];

		double X[3][3] = {
			{xyzCell[c0*3], xyzCell[c0*3+1], xyzCell[c0*3+2]},
			{xyzCell[c1*3], xyzCell[c1*3+1], xyzCell[c1*3+2]},
			{xyzCell[c2*3], xyzCell[c2*3+1], xyzCell[c2*3+2]},
		};

		double P[3] = {xVertex[i], yVertex[i], zVertex[i]};

		double lambda[3];
		barycentric_point2triangle(X[0], X[1], X[2], P, lambda);
		//[i * nVertLevels + curVertLevel];
		for (int curVertLevel=0;curVertLevel<nVertLevels; curVertLevel++){
			//        velocityXv[curVertLevel*nVertLevels+i] = lambda[0] * velocityX[c0*nVertLevels+curVertLevel] + lambda[1] * velocityX[c1*nVertLevels+curVertLevel] + lambda[2] * velocityX[c2*nVertLevels+curVertLevel];
			//        velocityYv[curVertLevel*nVertLevels+i] = lambda[0] * velocityY[c0*nVertLevels+curVertLevel] + lambda[1] * velocityY[c1*nVertLevels+curVertLevel] + lambda[2] * velocityY[c2*nVertLevels+curVertLevel];
			//        velocityZv[curVertLevel*nVertLevels+i] = lambda[0] * velocityZ[c0*nVertLevels+curVertLevel] + lambda[1] * velocityZ[c1*nVertLevels+curVertLevel] + lambda[2] * velocityZ[c2*nVertLevels+curVertLevel];

			velocityXv[i*nVertLevels+curVertLevel] = lambda[0] * velocityX[c0*nVertLevels+curVertLevel] + lambda[1] * velocityX[c1*nVertLevels+curVertLevel] + lambda[2] * velocityX[c2*nVertLevels+curVertLevel];
			velocityYv[i*nVertLevels+curVertLevel] = lambda[0] * velocityY[c0*nVertLevels+curVertLevel] + lambda[1] * velocityY[c1*nVertLevels+curVertLevel] + lambda[2] * velocityY[c2*nVertLevels+curVertLevel];
			velocityZv[i*nVertLevels+curVertLevel] = lambda[0] * velocityZ[c0*nVertLevels+curVertLevel] + lambda[1] * velocityZ[c1*nVertLevels+curVertLevel] + lambda[2] * velocityZ[c2*nVertLevels+curVertLevel];

			zTopVertex[nVertLevels*i + curVertLevel] = lambda[0] * zTop[c0*nVertLevels+curVertLevel] + lambda[1] * zTop[c1*nVertLevels+curVertLevel] + lambda[2] * zTop[c2*nVertLevels+curVertLevel];

		}

	}

	NC_SAFE_CALL( nc_close(ncid));


}

