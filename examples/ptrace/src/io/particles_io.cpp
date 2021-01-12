#include "particles_io.h"
#include <netcdf.h>
#include <pnetcdf.h>
#include "def.h"
#include "misc.h"


void load_MPAS_cells_coords(const std::string& filename, seed_data &data){

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

        size_t maxEdges = 6;

        NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

        NC_SAFE_CALL( nc_inq_dimid(ncid, "nCells", &dimid_cells) );
        NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_cells, &data.nCells) );


        NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
        NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
        NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );
        NC_SAFE_CALL( nc_inq_varid(ncid, "verticesOnCell", &varid_verticesOnCell));
        NC_SAFE_CALL( nc_inq_varid(ncid, "nEdgesOnCell", &varid_nEdgesOnCell));

        NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_vertices) );
        NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertices, &data.nVertices) );

        data.xVertex.resize(data.nVertices);
	    data.yVertex.resize(data.nVertices);
	    data.zVertex.resize(data.nVertices);

        size_t start_vertices[1] = {0}, size_vertices[1] = {data.nVertices};

        NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, &data.xVertex[0]) );
        NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, &data.yVertex[0]) );
        NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, &data.zVertex[0]) );

        
        size_t start_cell_maxedges[2] = {0,0}, size_cell_maxedges[2] = {data.nCells, maxEdges};
        data.verticesOnCell.resize(data.nCells * maxEdges);
        NC_SAFE_CALL( nc_get_vara_int(ncid, varid_verticesOnCell, start_cell_maxedges, size_cell_maxedges, &data.verticesOnCell[0]) );

        size_t start_cells[1] = {0}, size_cells[1] = {data.nCells};
        data.nEdgesOnCell.resize(data.nCells);
        NC_SAFE_CALL( nc_get_vara_int(ncid, varid_nEdgesOnCell, start_cells, size_cells, &data.nEdgesOnCell[0]) );


        NC_SAFE_CALL( nc_close(ncid));

}

void load_particle_data(const std::string& filename, seed_data &data){

    int ncid;

    int dimid_nStep, dimid_nParticle, dimid_nVertices, dimid_nVertexAllLayers;
    int varid_xPos, varid_yPos, varid_zPos, varid_tsg,
            varid_xVertex, varid_yVertex, varid_zVertex, varid_zTopVertex,
            varid_velocityVx, varid_velocityVy, varid_velocityVz; // trace_sizes_global
    int varid_zLevelParticle;
    size_t nParticles, nSteps;



    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "Time", &dimid_nStep) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "nParticles", &dimid_nParticle) );
//    NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_nVertices) );
//    NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertexAllLayers", &dimid_nVertexAllLayers) );

    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nParticle, &nParticles) );
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nStep, &nSteps) );
//    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nVertices, &nVertices) );
//    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nVertexAllLayers, &nVertexAllLayers) );

    fprintf(stderr, "particles.cpp load_particle_data %ld %ld\n", nParticles, nSteps);

    NC_SAFE_CALL( nc_inq_varid(ncid, "xParticle", &varid_xPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "yParticle", &varid_yPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "zParticle", &varid_zPos) );




    const size_t start_pos[2] = {0, 0};
    const size_t count_pos[2] = {nSteps, nParticles};

    data.xPos.resize(nSteps * nParticles);
    data.yPos.resize(nSteps * nParticles);
    data.zPos.resize(nSteps * nParticles);
    // zLevelParticle.resize(nSteps * nParticles);
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xPos, start_pos, count_pos, &data.xPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yPos, start_pos, count_pos, &data.yPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zPos, start_pos, count_pos, &data.zPos[0]) );
    // NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zLevelParticle, start_pos, count_pos, &zLevelParticle[0]) );



    NC_SAFE_CALL( nc_close(ncid));


}