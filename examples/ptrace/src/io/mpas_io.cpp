#include "mpas_io.h"
#include "misc.h"
#include <netcdf.h>
#include <pnetcdf.h>
#include "def.h"

mpas_io::mpas_io(){

		// velocityXv.resize(2);
		// velocityYv.resize(2);
		// velocityZv.resize(2);
		// zTop.resize(2);
		// zMid.resize(2);
		// vertVelocityTop.resize(2);


		// TODO: replace hard code values, also by local sizes
		nCells = 7234;
		nEdges =  22736;
		nVertices = 15459;
		nVertLevels = 100;
		maxEdges = 6;
		nVertLevelsP1 = 101;
		radius = 6371229.;
		cradius = radius;

		velocitiesX.resize(2);
		velocitiesY.resize(2);
		velocitiesZ.resize(2);


}

void mpas_io::move_data_to_regular_position(int gid){

	// create buffers
	std::vector<double> tmp_velocityX(nCells * nVertLevels);
	std::vector<double> tmp_velocityY(nCells * nVertLevels);
	std::vector<double> tmp_velocityZ(nCells * nVertLevels);
	std::vector<double> tmp_vertVelocityTop(nCells * nVertLevelsP1);
	std::vector<double> tmp_zTop(nCells * nVertLevels);
	std::vector<double> tmp_zMid(nCells * nVertLevels);
	
	// dprint("indexToCellid size %ld | nC %ld, nVL %ld, nVLP1 %ld", indexToCellID.size(), nCells, nVertLevels, nVertLevelsP1);

	// dprint("velocitiesY[1] %ld, nCells %ld nVertLevels %ld indexToCellID %ld", velocitiesY[1].size(), nCells, nVertLevels, indexToCellID.size());

	// iterate through each cell and copy data to tmp
	int ctr = 0;
	int flag = 0;
	for (auto gcId: indexToCellID){

		// local_gcIds_init.push_back(id);

		int idx = gcId - 1;

		

		// if (id == 2630)	{
		// 	dprint("found id %d , ctr %d", id, ctr);
		// 	flag =1;
		// 	for (size_t i=ctr*nVertLevels; i<ctr*nVertLevels+nVertLevels; i++)
		// 		fprintf(stderr, "%f ", velocityX[i]);
		// }

		

		

		std::copy(velocitiesX[1].begin()+ ctr*nVertLevels, velocitiesX[1].begin()+ ctr*nVertLevels + nVertLevels, tmp_velocityX.begin() + idx * nVertLevels);
		std::copy(velocitiesY[1].begin()+ ctr*nVertLevels, velocitiesY[1].begin()+ ctr*nVertLevels + nVertLevels, tmp_velocityY.begin() + idx * nVertLevels);
		std::copy(velocitiesZ[1].begin()+ ctr*nVertLevels, velocitiesZ[1].begin()+ ctr*nVertLevels + nVertLevels, tmp_velocityZ.begin() + idx * nVertLevels);
		std::copy(vertVelocityTop.begin()+ ctr*nVertLevelsP1, vertVelocityTop.begin()+ ctr*nVertLevelsP1 + nVertLevelsP1, tmp_vertVelocityTop.begin()+idx*nVertLevelsP1);
		std::copy(zTop.begin()+ ctr*nVertLevels, zTop.begin()+ ctr*nVertLevels + nVertLevels, tmp_zTop.begin() + idx * nVertLevels);
		std::copy(zMid.begin()+ ctr*nVertLevels, zMid.begin()+ ctr*nVertLevels + nVertLevels, tmp_zMid.begin() + idx * nVertLevels);

		ctr++;
	}

	// dprint("later indexToCellid size %ld", indexToCellID.size());

	// move the tmp to mpas1
	velocitiesX[1] = std::move(tmp_velocityX);
	velocitiesY[1] = std::move(tmp_velocityY);
	velocitiesZ[1] = std::move(tmp_velocityZ);
	vertVelocityTop = std::move(tmp_vertVelocityTop);
	zTop = std::move (tmp_zTop);
	zMid = std::move(tmp_zMid);

	if (velocitiesX[0].size()==0){
		velocitiesX[0].resize(velocitiesX[1].size());
		velocitiesY[0].resize(velocitiesY[1].size());
		velocitiesZ[0].resize(velocitiesZ[1].size());
	}

	// ctr = 2630 - 1;
	// if (gid==0)
	// 	for (size_t i=ctr*nVertLevels; i<ctr*nVertLevels + nVertLevels; i++)
	// 		fprintf(stderr, "%f ", velocitiesX[1][i]);
	

		
}

void mpas_io::update_data(int data_id, int frame_no, std::vector<int> &data_int, std::vector<double> &data_dbl){


	

	switch(data_id){

				case 0: //fprintf(stderr, "Recv xCell %d,\n", data_id);
						// xCell = std::move(data_dbl);
						break;

				case 1: //fprintf(stderr, "Recv yCell %d,\n", data_id);
						// yCell = std::move(data_dbl);
						break;

				case 2:	
						//printf(stderr, "Recv zCell %d,\n", data_id);
						// zCell = std::move(data_dbl);
						break;

				case 3:	//fprintf(stderr, "Recv xVertex %d,\n", data_id);
						//xVertex = std::move(data_dbl);
						break;
						
				case 4:	//fprintf(stderr, "Recv yVertex %d,\n", data_id);
						//yVertex = std::move(data_dbl);
						break;

				case 5:	//fprintf(stderr, "Recv zVertex %d,\n", data_id);
						//zVertex = std::move(data_dbl);
						break;

				case 6:	//fprintf(stderr, "Recv indexToVertexID %d,\n", data_id);
						//indexToVertexID = std::move(data_int);

						// for (int i=0; i<indexToVertexID.size(); i++) { 
						// 		vertexIndex[indexToVertexID[i]] = i;
						// }
						
						break;

				case 7:	//fprintf(stderr, "Recv indexToCellID %d,\n", data_id);
						indexToCellID = std::move(data_int);
						
						// for (int i=0; i<indexToCellID.size(); i++) { 
						// 		cellIndex[indexToCellID[i]] = i;
						// }
						break;

				case 8:	//fprintf(stderr, "Recv verticesOnEdge %d,\n", data_id);
						//verticesOnEdge = std::move(data_int);
						break;

				case 9:	//fprintf(stderr, "Recv cellsOnVertex %d,\n", data_id);
						//cellsOnVertex = std::move(data_int);
						break;

				case 10:{	//fprintf(stderr, "Recv verticesOnCell %d,\n", data_id);
						//verticesOnCell = std::move(data_int);
						// size_t max=0;
						// 	for (size_t i=0; i<verticesOnCell.size(); i++){
						// 		if (max<verticesOnCell[i])
						// 			max = verticesOnCell[i];
						// 	}
							// dprint("verticesOnCell %ld", verticesOnCell.size());
						}
						break;
						

				case 11: 	//fprintf(stderr, "Recv velocityXv %d, fno %d\n", data_id, frame_no);
						velocitiesX[0] = std::move(velocitiesX[1]);
						velocitiesX[1] = std::move(data_dbl);
						
						// dprint("first velocityXv %f", velocityXv[frame_no%2][0]);
						break;

				case 12:	//fprintf(stderr, "Recv velocityYv %d, fno %d\n", data_id, frame_no);
						velocitiesY[0] = std::move(velocitiesY[1]);
						velocitiesY[1] = std::move(data_dbl);
						
						break;

				case 13:	//fprintf(stderr, "Recv velocityZv %d,\n", data_id);
						velocitiesZ[0] = std::move(velocitiesZ[1]);
						velocitiesZ[1] = std::move(data_dbl);
						
						break;

				case 14://	fprintf(stderr, "Recv nEdgesOnCell %d,\n", data_id);
						// nEdgesOnCell = std::move(data_int);
						break;

				case 15://	fprintf(stderr, "Recv maxLevelCell %d,\n", data_id);
						// maxLevelCell = std::move(data_int);
						break;

				case 16:	//fprintf(stderr, "Recv boundaryVertex %d,\n", data_id);
						// boundaryVertex = std::move(data_int);
						break;

				case 17:	//fprintf(stderr, "Recv vertVelocityTop %d,\n", data_id);
						vertVelocityTop = std::move(data_dbl);
						// dprint("vertVelocityTop.. %f %f", vertVelocityTop[frame_no%2][0], vertVelocityTop[frame_no%2][1]);
						break;

				case 18:{//fprintf(stderr, "Recv cellsOnCell %d,\n", data_id);
						 //	cellsOnCell = std::move(data_int);
						} 
						
						break;

				case 19:	//fprintf(stderr, "Recv zTop %d,\n", data_id);
						zTop = std::move(data_dbl);
						break;

				case 20:	//fprintf(stderr, "Recv zMid %d,\n", data_id);
						zMid = std::move(data_dbl);
						break;

				default:fprintf(stderr, "Exiting! %d,\n", data_id); 
						exit(0);
						break;
				}

}

void mpas_io::create_cells_kdtree(){


		int ncid;
		int dimid_cells;
		int varid_xCell, varid_yCell, varid_zCell;
		size_t nCells;
		

		// std::vector<int> glCellIdx;

		// output2 to deal with online updation of output.nc by simulation
		NC_SAFE_CALL( nc_open("output2.nc", NC_NOWRITE, &ncid) );

		NC_SAFE_CALL( nc_inq_dimid(ncid, "nCells", &dimid_cells) );
		NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_cells, &nCells) );

		

		NC_SAFE_CALL( nc_inq_varid(ncid, "xCell", &varid_xCell) );
		NC_SAFE_CALL( nc_inq_varid(ncid, "yCell", &varid_yCell) );
		NC_SAFE_CALL( nc_inq_varid(ncid, "zCell", &varid_zCell) );

		dprint("nCells %ld", nCells);

		std::vector<double> xCells, yCells, zCells;
		xCells.resize(nCells);
		yCells.resize(nCells);
		zCells.resize(nCells);

		const size_t start_cells[1] = {0}, size_cells[1] = {nCells};
		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &xCells[0]) );
		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &yCells[0]) );
		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &zCells[0]) );

		
		// dprint("varids %d %d %d,  %f %f %f", varid_xCell, varid_yCell, varid_zCell, xCells[1], yCells[1], zCells[1]);

		// dprint("here %ld", xCells.size());

		NC_SAFE_CALL( nc_close(ncid));

		C_local.resize(3, xCells.size());
		for (int i=0;i<xCells.size();i++){
			C_local(0,i) = xCells[i];
			C_local(1,i) = yCells[i];
			C_local(2,i) = zCells[i];

      // if (i==23 || i==2298){
      //   dprint("i %d, %f %f %f ", i, xyzCell[i*3], xyzCell[i*3+1], xyzCell[i*3+2]);
      // }
		}

		nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);
		


}


void mpas_io::loadMeshFromNetCDF_CANGA(diy::mpi::communicator& world, const std::string& filename, long long int time_id, const std::string &fname_graph){


	int ncid;
	int dimid_cells, dimid_edges, dimid_vertices, dimid_vertLevels, dimid_vertLevelsP1;
	int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex,
	    varid_latCell, varid_lonCell, varid_xCell, varid_yCell, varid_zCell,
	    varid_verticesOnEdge, varid_cellsOnVertex,
	    varid_indexToVertexID, varid_indexToCellID,
	    varid_velocityX, varid_velocityY, varid_velocityZ, 
		varid_uVertexVelocity, varid_vVertexVelocity, varid_wVertexVelocity, 
		varid_cellsOnCell, varid_zMid, varid_vertVelocityTop, varid_nEdgesOnCell, 
		varid_maxLevelCell, varid_verticesOnCell, varid_boundaryVertex;
		
	int varid_zTop;

	PNC_SAFE_CALL( ncmpi_open(world, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid) );

	PNC_SAFE_CALL( ncmpi_inq_dimid(ncid, "nCells", &dimid_cells) );
	PNC_SAFE_CALL( ncmpi_inq_dimid(ncid, "nEdges", &dimid_edges) );
	PNC_SAFE_CALL( ncmpi_inq_dimid(ncid, "nVertices", &dimid_vertices) );
	PNC_SAFE_CALL( ncmpi_inq_dimid(ncid, "nVertLevels", &dimid_vertLevels) );
	PNC_SAFE_CALL( ncmpi_inq_dimid(ncid, "nVertLevelsP1", &dimid_vertLevelsP1) );

	PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimid_cells, &nCells) );
	PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimid_edges, &nEdges) );
	PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimid_vertices, &nVertices) );
	PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimid_vertLevels, &nVertLevels) );
	PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimid_vertLevelsP1, &nVertLevelsP1) );

	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "indexToVertexID", &varid_indexToVertexID) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "indexToCellID", &varid_indexToCellID) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "latCell", &varid_latCell) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "lonCell", &varid_lonCell) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "xCell", &varid_xCell) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "yCell", &varid_yCell) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "zCell", &varid_zCell) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "latVertex", &varid_latVertex) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "lonVertex", &varid_lonVertex) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "xVertex", &varid_xVertex) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "yVertex", &varid_yVertex) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "zVertex", &varid_zVertex) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "verticesOnEdge", &varid_verticesOnEdge) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "cellsOnVertex", &varid_cellsOnVertex) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "velocityX", &varid_velocityX) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "velocityY", &varid_velocityY) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "velocityZ", &varid_velocityZ) );
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "zTop", &varid_zTop));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "boundaryVertex", &varid_boundaryVertex));

	// PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "uVertexVelocity", &varid_uVertexVelocity));
	// PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "vVertexVelocity", &varid_vVertexVelocity));
	// PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "wVertexVelocity", &varid_wVertexVelocity));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "cellsOnCell", &varid_cellsOnCell));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "zMid", &varid_zMid));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "vertVelocityTop", &varid_vertVelocityTop));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "nEdgesOnCell", &varid_nEdgesOnCell));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "maxLevelCell", &varid_maxLevelCell));
	PNC_SAFE_CALL( ncmpi_inq_varid(ncid, "verticesOnCell", &varid_verticesOnCell));


	const MPI_Offset start_cells[1] = {0}, size_cells[1] = {nCells};

	// indexToCellID.resize(nCells);
	// PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_indexToCellID, start_cells, size_cells, &indexToCellID[0]) );
	// for (int i=0; i<nCells; i++) {
	// 	cellIndex[indexToCellID[i]] = i;
	// 	// fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
	// }

	std::vector<double> coord_cells;
	coord_cells.resize(nCells);
	xyzCell.resize(nCells*3);

	xCell.resize(nCells);
	yCell.resize(nCells);
	zCell.resize(nCells);

	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_xCell, start_cells, size_cells, &xCell[0]) );
	for (int i=0; i<nCells; i++) xyzCell[i*3] = xCell[i];
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_yCell, start_cells, size_cells, &yCell[0]) );
	for (int i=0; i<nCells; i++) xyzCell[i*3+1] = yCell[i];
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_zCell, start_cells, size_cells, &zCell[0]) );
	for (int i=0; i<nCells; i++) xyzCell[i*3+2] = zCell[i];

	maxLevelCell.resize(nCells);
	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_maxLevelCell, start_cells, size_cells, &maxLevelCell[0]) );	

	const MPI_Offset start_vertices[1] = {0}, size_vertices[1] = {nVertices};
	latVertex.resize(nVertices);
	lonVertex.resize(nVertices);
	xVertex.resize(nVertices);
	yVertex.resize(nVertices);
	zVertex.resize(nVertices);
	indexToVertexID.resize(nVertices);

	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_indexToVertexID, start_vertices, size_vertices, &indexToVertexID[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_latVertex, start_vertices, size_vertices, &latVertex[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_lonVertex, start_vertices, size_vertices, &lonVertex[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_xVertex, start_vertices, size_vertices, &xVertex[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_yVertex, start_vertices, size_vertices, &yVertex[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_zVertex, start_vertices, size_vertices, &zVertex[0]) );

	// for (int i=0; i<nVertices; i++) {
	// 	vertexIndex[indexToVertexID[i]] = i;
	// 	// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	// }

	const MPI_Offset start_edges2[2] = {0, 0}, size_edges2[2] = {nEdges, 2};
	verticesOnEdge.resize(nEdges*2);

	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_verticesOnEdge, start_edges2, size_edges2, &verticesOnEdge[0]) );

	// for (int i=0; i<nEdges; i++)
	//   fprintf(stderr, "%d, %d\n", verticesOnEdge[i*2], verticesOnEdge[i*2+1]);


	const MPI_Offset start_vertex_cell[2] = {0, 0}, size_vertex_cell[2] = {nVertices, 3};
	cellsOnVertex.resize(nVertices*3);

	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_cellsOnVertex, start_vertex_cell, size_vertex_cell, &cellsOnVertex[0]) );

	// getting velocity at one time at cell centers
	const MPI_Offset start_time_cell_level[3] = {time_id, 0, 0}, size_time_cell_level[3] = {1, nCells, nVertLevels};
	velocityX.resize(nCells*nVertLevels);
	velocityY.resize(nCells*nVertLevels);
	velocityZ.resize(nCells*nVertLevels);

	// dprint("time_id %ld, nCells %ld nVertLevels %ld", time_id, nCells, nVertLevels);



	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_velocityX, start_time_cell_level, size_time_cell_level, &velocityX[0]) );
	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_velocityY, start_time_cell_level, size_time_cell_level, &velocityY[0]) );
	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_velocityZ, start_time_cell_level, size_time_cell_level, &velocityZ[0]) );

	zTop.resize(nCells*nVertLevels);
	zMid.resize(nCells*nVertLevels);
	zTopVertex.resize(nVertices*nVertLevels);
	zTopVertexNorm.resize(nVertices*nVertLevels);

	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_zTop, start_time_cell_level, size_time_cell_level, &zTop[0]) );
	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_zMid, start_time_cell_level, size_time_cell_level, &zMid[0]) );
	vertVelocityTop.resize(nCells * nVertLevelsP1);
	const MPI_Offset start_time_cell_levelP1[3] = {time_id, 0, 0}, size_time_cell_levelP1[3] = {1, nCells, nVertLevelsP1};
	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_vertVelocityTop, start_time_cell_level, size_time_cell_level, &vertVelocityTop[0]) );


	// getting velocity at one time at vertex centers
	const MPI_Offset start_time_vertex_level[3] = {time_id, 0, 0}, size_time_vertex_level[3] = {1, nVertices, nVertLevels};

	// derive velocity on verticies (on only one time level here?)
	uVertexVelocity.resize(nVertices*nVertLevels);
	vVertexVelocity.resize(nVertices*nVertLevels);
	wVertexVelocity.resize(nVertices*nVertLevels);

	uVertexVelocities.resize(2); 
	uVertexVelocities[0].resize(nVertices*nVertLevels);	
	uVertexVelocities[1].resize(nVertices*nVertLevels);	

	vVertexVelocities.resize(2);
	vVertexVelocities[0].resize(nVertices*nVertLevels);	
	vVertexVelocities[1].resize(nVertices*nVertLevels);	

	wVertexVelocities.resize(2);
	wVertexVelocities[0].resize(nVertices*nVertLevels);	
	wVertexVelocities[1].resize(nVertices*nVertLevels);	

	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_uVertexVelocity, start_time_vertex_level, size_time_vertex_level, &uVertexVelocity[0]) );
	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_vVertexVelocity, start_time_vertex_level, size_time_vertex_level, &vVertexVelocity[0]) );
	// PNC_SAFE_CALL( ncmpi_get_vara_double_all(ncid, varid_wVertexVelocity, start_time_vertex_level, size_time_vertex_level, &wVertexVelocity[0]) );
	cellsOnCell.resize(nCells * maxEdges);
	nEdgesOnCell.resize(nCells);
	verticesOnCell.resize(nCells * maxEdges);
	const MPI_Offset start_cell_maxedges[2] = {0,0}, size_cell_maxedges[2] = {nCells, maxEdges};
	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_cellsOnCell, start_cell_maxedges, size_cell_maxedges, &cellsOnCell[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_nEdgesOnCell, start_cells, size_cells, &nEdgesOnCell[0]) );
	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_verticesOnCell, start_cell_maxedges, size_cell_maxedges, &verticesOnCell[0]) );
	
	boundaryVertex.resize(nVertices * nVertLevels);
	const MPI_Offset start_nvert_nvertlevel[2] = {0, 0}, size_nvert_nvertlevel[2] = {nVertices, nVertLevels};
	PNC_SAFE_CALL( ncmpi_get_vara_int_all(ncid, varid_boundaryVertex, start_nvert_nvertlevel, size_nvert_nvertlevel, &boundaryVertex[0]) );


	PNC_SAFE_CALL( ncmpi_close(ncid));

	cell_nbrs = read_csv(fname_graph, ' ');

	// int vert = 43382;
	// dprint("vertOnCell %d %d %d %d %d %d", verticesOnCell[vert*6], verticesOnCell[vert*6+1], verticesOnCell[vert*6+2], verticesOnCell[vert*6+3], verticesOnCell[vert*6+4], verticesOnCell[vert*6+5]);

	// for (size_t i=0; i<verticesOnCell.size(); i++){
	// 	if (verticesOnCell[i]==429305 && world.rank()==0)
	// 		dprint("cellID %ld %ld", i, i/6);
	// }

	// vertOnCell 434468 434488 429306 430222 430708 429126
	// int vid = 434468-1; 
	// dprint("uVertVel %f %f %f", uVertexVelocity[vid*nVertLevels], uVertexVelocity[vid*nVertLevels+1], uVertexVelocity[vid*nVertLevels+2]);
	// int cell = 43382;
	// dprint("cell coords %f %f %f", xCell[cell], yCell[cell], zCell[cell]);

	// dprint("nVertL %ld, nCells %ld, nVertLP1 %ld", nVertLevels, nCells, nVertLevelsP1);


}