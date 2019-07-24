
#include "block.h"
#include "misc.h"
#include <netcdf.h>
#include "def.h"
#include <Eigen/Dense>
#include "nabo/nabo.h"
#include <map>
#include <utility>




std::vector<Halo> block::get_halo_info(){

	std::map<int, std::pair<std::vector<int>, std::vector<int>>> cells_verts;

	// iterate over cells
	for(size_t i=0; i<indexToCellID.size(); i++){

		int cellID = indexToCellID[i];
		// iterate over cell neighbors
		for (int j = i*maxEdges; j<(i+1)*maxEdges; j++){

			int nbr_cellID = cellsOnCell[j];
			// make sure nbr_cell is valid
			if (nbr_cellID != 0){
				int nbr_gid = gcIdxToGid[nbr_cellID-1];
				// if neighbor gid != current gid add then add cell id request to nbr_gid 
					if (nbr_gid != gid){
							cells_verts[nbr_gid].first.push_back(nbr_cellID);
					}
			}

		}

		// iterate over vertices on cell
		for (size_t j=i*maxEdges; j<(i+1)*maxEdges; j++){

			int nbr_vertexID = verticesOnCell[j];

			// make sure nbr_vertex is valid
			if (nbr_vertexID !=0 ){
				int nbr_gid = gVIdxToGid[nbr_vertexID-1];

				// if neighbor gid != current gid add then add vertex id request to nbr_gid 
				if (nbr_gid != gid){
					cells_verts[nbr_gid].second.push_back(nbr_vertexID);
				}

			}

			

		}
	}


	// iterate over map and populate halos

	std::vector<Halo> halos;

	std::for_each(cells_verts.begin(), cells_verts.end(),
		[&](std::pair<int, std::pair<std::vector<int>, std::vector<int>>> element){

			Halo h;
			h.src_gid = gid;
			h.dest_gid = element.first;
			h.glCellIDs = std::move(element.second.first);
			h.glVertexIDs = std::move(element.second.second);

			halos.push_back(h);
		});

	return halos;
}


void block::process_halo_req(Halo &h, int framenum){

	// store the shell to populate dynamic halo containers later for velocityXv,...
	halo_info.push_back(h);

	// populate cell info for requested cellids
	for (size_t i=0; i<h.glCellIDs.size(); i++){
		h.xCell.push_back(xCell[cellIndex[h.glCellIDs[i]]]);
		h.yCell.push_back(yCell[cellIndex[h.glCellIDs[i]]]);
		h.zCell.push_back(zCell[cellIndex[h.glCellIDs[i]]]);
	}

	// populate vert info for requested vertex ids

	for (size_t i=0; i<h.glVertexIDs.size(); i++){

		h.xVertex.push_back(xVertex[vertexIndex[h.glVertexIDs[i]]]);
		h.yVertex.push_back(yVertex[vertexIndex[h.glVertexIDs[i]]]);
		h.zVertex.push_back(zVertex[vertexIndex[h.glVertexIDs[i]]]);
		for(size_t j=0; j<nVertLevels; j++){
			h.velocityXv.push_back(velocityXv[framenum%2][j+nVertLevels*vertexIndex[h.glVertexIDs[i]]]);
			h.velocityYv.push_back(velocityYv[framenum%2][j+nVertLevels*vertexIndex[h.glVertexIDs[i]]]);
			h.velocityZv.push_back(velocityZv[framenum%2][j+nVertLevels*vertexIndex[h.glVertexIDs[i]]]);
			h.vertVelocityTop.push_back(vertVelocityTop[framenum%2][j+nVertLevels*vertexIndex[h.glVertexIDs[i]]]);
			h.zMid.push_back(zMid[framenum%2][j+nVertLevels*vertexIndex[h.glVertexIDs[i]]]);
			h.zTop.push_back(zTop[framenum%2][j+nVertLevels*vertexIndex[h.glVertexIDs[i]]]);

		}
	}

}


void block::update_halo_info(Halo &h, int framenum){

	

	// iterate through each incoming cell and append info to local arrays

	size_t nCells_local = xCell.size();
	size_t nCells_all = nCells_local+h.glCellIDs.size();

	xCell.resize(nCells_all);
	yCell.resize(nCells_all);
	zCell.resize(nCells_all);

	for (size_t i=0; i<h.glCellIDs.size(); i++){

		int cellID = h.glCellIDs[i];

		cellIndex[cellID] = indexToCellID.size();
		indexToCellID.push_back(cellID);
		
		xCell[nCells_local+i] = h.xCell[i];
		yCell[nCells_local+i] = h.yCell[i];
		zCell[nCells_local+i] = h.zCell[i];
	}




	// iterate through each incoming vertex and append info to local arrays

	size_t nVertices_local = xVertex.size();
	size_t nVertices_all = nVertices_local+h.glVertexIDs.size();

	xVertex.resize(nVertices_all);
	yVertex.resize(nVertices_all);
	zVertex.resize(nVertices_all);
	velocityXv[framenum%2].resize(nVertices_all*nVertLevels);
	velocityYv[framenum%2].resize(nVertices_all*nVertLevels);
	velocityZv[framenum%2].resize(nVertices_all*nVertLevels);
	zMid[framenum%2].resize(nVertices_all*nVertLevels);
	zTop[framenum%2].resize(nVertices_all*nVertLevels);

	for (size_t i=0; i<h.glVertexIDs.size(); i++){

		int vertID = h.glVertexIDs[i];
		vertexIndex[vertID] = indexToVertexID.size();
		indexToVertexID.push_back(vertID);

		xVertex[nVertices_local+i] = h.xVertex[i];
		yVertex[nVertices_local+i] = h.yVertex[i];
		zVertex[nVertices_local+i] = h.zVertex[i];

		for(size_t j=0; j<nVertLevels; j++){
			velocityXv[framenum%2][j+nVertLevels*(nVertices_local+i)] = h.velocityXv[j+nVertLevels*i];
			velocityYv[framenum%2][j+nVertLevels*(nVertices_local+i)] = h.velocityYv[j+nVertLevels*i];
			velocityZv[framenum%2][j+nVertLevels*(nVertices_local+i)] = h.velocityZv[j+nVertLevels*i];
			zMid[framenum%2][j+nVertLevels*(nVertices_local+i)] = h.zTop[j+nVertLevels*i];
			zTop[framenum%2][j+nVertLevels*(nVertices_local+i)] = h.zTop[j+nVertLevels*i];

		}

	}


}

void block::process_halo_dynamic(Halo &h, int framenum){


}


void block::update_halo_dynamic(Halo &h, int framenum){


	
}

void block::create_links_mpas(const std::string &fname_graphinfo, std::set<int> &links, diy::mpi::communicator &world){

	std::vector<std::vector<int>> partn_ids = read_csv(fname_graphinfo.c_str());

	for (size_t i=0; i<partn_ids.size(); i++)
		gcIdxToGid.push_back(partn_ids[i][0]);

	for (size_t i=0; i<partn_ids.size(); i++){
		if (gid==partn_ids[i][0]){
			for (int j = cellIndex[i+1]*maxEdges; j<(cellIndex[i+1]+1)*maxEdges; j++){
				if (cellsOnCell[j]!=0){
					if (partn_ids[cellsOnCell[j]-1][0] != gid){
						links.insert(partn_ids[cellsOnCell[j]-1][0]);
					}
				}
			}
		}
	}		

	// dprint("links size %ld", links.size());

	/* now populate gVIdxtoGid */

	// get nVertices_global
	diy::mpi::all_reduce(world, indexToVertexID.size(), nVertices_global, std::plus<size_t>());
	gVIdxToGid.resize(nVertices_global);

	// all_gather the gids
	std::vector<int> all_gids;
	diy::mpi::all_gather(world, gid, all_gids);


	// all gather indexToVertexID
	std::vector<std::vector<int>> all_indexToVertexID;
	diy::mpi::all_gather(world, indexToVertexID, all_indexToVertexID);
	for (size_t i=0; i<all_indexToVertexID.size(); i++){
		for (size_t j=0; j<all_indexToVertexID[i].size(); j++){

			gVIdxToGid[all_indexToVertexID[i][j]-1] = all_gids[i];
		}
	}


}


// create a new particle file which includes a glCellIdx (global cell index) field
void block::generate_new_particle_file(){

	if (gid==0){

		// read output.nc & create kd-tree

		int ncid;
		int dimid_cells;
		int varid_xCell, varid_yCell, varid_zCell;
		size_t nCells;
		

		std::vector<int> glCellIdx;

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

		NC_SAFE_CALL( nc_close(ncid));


		
		/* KDTree */
		Eigen::MatrixXd C_local;
		Nabo::NNSearchD*  nns_cells;

		C_local.resize(3, nCells);
		for (int i=0;i<nCells;i++){
			C_local(0,i) = xCells[i];
			C_local(1,i) = yCells[i];
			C_local(2,i) = zCells[i];
		}

		nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);


	// read particles.nc

		int dimid_particles, dimid_time;
		int varid_xParticle, varid_yParticle, varid_zParticle, varid_zLevelParticle;
		std::vector<double> xParticle, yParticle, zParticle, zLevelParticle;
		size_t nParticles;
		size_t nTime;

		NC_SAFE_CALL(nc_open("particles.nc", NC_CLOBBER, &ncid));

		NC_SAFE_CALL( nc_inq_dimid(ncid, "nParticles", &dimid_particles) );
		NC_SAFE_CALL( nc_inq_dimid(ncid, "Time", &dimid_time) );

		NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_particles, &nParticles));
		NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_time, &nTime));

		NC_SAFE_CALL( nc_inq_varid(ncid, "xParticle", &varid_xParticle));
		NC_SAFE_CALL( nc_inq_varid(ncid, "yParticle", &varid_yParticle));
		NC_SAFE_CALL( nc_inq_varid(ncid, "zParticle", &varid_zParticle));
		NC_SAFE_CALL( nc_inq_varid(ncid, "zLevelParticle", &varid_zLevelParticle));

		xParticle.resize(nTime*nParticles);
		yParticle.resize(nTime*nParticles);
		zParticle.resize(nTime*nParticles);
		zLevelParticle.resize(nTime*nParticles);
		glCellIdx.resize(nParticles);
		const size_t start_t_p[2] = {0, 0}, size_t_p[2] = {nTime, nParticles};

		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xParticle, start_t_p, size_t_p, &xParticle[0]));
		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yParticle, start_t_p, size_t_p, &yParticle[0]));
		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zParticle, start_t_p, size_t_p, &zParticle[0]));
		NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zLevelParticle, start_t_p, size_t_p, &zLevelParticle[0]));

		NC_SAFE_CALL(nc_close(ncid));


	// parce, filter and pushback particles



		int cellID;
		for (size_t i = 0; i<xParticle.size(); i++){
			Eigen::VectorXd q(3);
			q<<xParticle[i], yParticle[i], zParticle[i];
        // get nearest cell neighbor ids using mpas_c
			Eigen::VectorXi nearest_cell_idx(1);
			Eigen::VectorXd dists2_cell(1);
			nns_cells->knn(q, nearest_cell_idx,dists2_cell, 1,0, Nabo::NNSearchF::SORT_RESULTS| Nabo::NNSearchF::ALLOW_SELF_MATCH);
        cellID = nearest_cell_idx[0] ; // cell global idx
        glCellIdx[i] = cellID;
        dprint("i %ld cellID %d, %f %f %f", i, cellID, xCells[i], yCells[i], zCells[i]);
    }



    /* Write particle2.nc file */

    int varid_glCellIdx;

    NC_SAFE_CALL(nc_create("particles2.nc", NC_CLOBBER, &ncid));

    NC_SAFE_CALL(nc_def_dim(ncid, "nParticles", nParticles, &dimid_particles));
    NC_SAFE_CALL(nc_def_dim(ncid, "Time", NC_UNLIMITED, &dimid_time));


    /* Define variables */
    int ndims_particle = 2;
    int dimids_particle[] = {dimid_time, dimid_particles};
    NC_SAFE_CALL(nc_def_var(ncid, "xParticle", NC_DOUBLE, ndims_particle, dimids_particle, &varid_xParticle));
    NC_SAFE_CALL(nc_def_var(ncid, "yParticle", NC_DOUBLE, ndims_particle, dimids_particle, &varid_yParticle));
    NC_SAFE_CALL(nc_def_var(ncid, "zParticle", NC_DOUBLE, ndims_particle, dimids_particle, &varid_zParticle));
    NC_SAFE_CALL(nc_def_var(ncid, "zLevelParticle", NC_DOUBLE, ndims_particle, dimids_particle, &varid_zLevelParticle));

    int dimids_glCellIdx[] = {dimid_particles};
    NC_SAFE_CALL(nc_def_var(ncid, "glCellIdx", NC_INT, 1, dimids_glCellIdx, &varid_glCellIdx));

    NC_SAFE_CALL(nc_enddef(ncid));

    size_t start[2], count[2];
    start[0] = 0; start[1] = 0;  count[0]=nTime; count[1] = nParticles;
    NC_SAFE_CALL(nc_put_vara_double(ncid, varid_xParticle, start, count, &xParticle[0]));
    NC_SAFE_CALL(nc_put_vara_double(ncid, varid_yParticle, start, count, &yParticle[0]));
    NC_SAFE_CALL(nc_put_vara_double(ncid, varid_zParticle, start, count, &zParticle[0]));
    NC_SAFE_CALL(nc_put_vara_double(ncid, varid_zLevelParticle, start, count, &zLevelParticle[0]));

    size_t start_glCellIdx[1]={0}, count_glCellIdx[1]={nParticles};
    NC_SAFE_CALL(nc_put_vara_int(ncid, varid_glCellIdx, start_glCellIdx, count_glCellIdx, &glCellIdx[0]));
    

    NC_SAFE_CALL(nc_close(ncid));


} // if (gid==0)



}

void block::init_seeds_mpas(std::string &fname_particles){

	
	// read particles2.nc
	int ncid;
	int dimid_particles, dimid_time;
	int varid_xParticle, varid_yParticle, varid_zParticle, varid_zLevelParticle, varid_glCellIdx;
	std::vector<double> xParticle, yParticle, zParticle, zLevelParticle;
	std::vector<int> glCellIdx;

	size_t nParticles;
	size_t nTime;

	NC_SAFE_CALL(nc_open(fname_particles.c_str(), NC_CLOBBER, &ncid));

	NC_SAFE_CALL( nc_inq_dimid(ncid, "nParticles", &dimid_particles) );
	NC_SAFE_CALL( nc_inq_dimid(ncid, "Time", &dimid_time) );

	NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_particles, &nParticles));
	NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_time, &nTime));

	NC_SAFE_CALL( nc_inq_varid(ncid, "xParticle", &varid_xParticle));
	NC_SAFE_CALL( nc_inq_varid(ncid, "yParticle", &varid_yParticle));
	NC_SAFE_CALL( nc_inq_varid(ncid, "zParticle", &varid_zParticle));
	NC_SAFE_CALL( nc_inq_varid(ncid, "zLevelParticle", &varid_zLevelParticle));
	NC_SAFE_CALL( nc_inq_varid(ncid, "glCellIdx", &varid_glCellIdx));

	xParticle.resize(nTime*nParticles);
	yParticle.resize(nTime*nParticles);
	zParticle.resize(nTime*nParticles);
	zLevelParticle.resize(nTime*nParticles);
	glCellIdx.resize(nParticles);
	const size_t start_t_p[2] = {0, 0}, size_t_p[2] = {nTime, nParticles};

	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xParticle, start_t_p, size_t_p, &xParticle[0]));
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yParticle, start_t_p, size_t_p, &yParticle[0]));
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zParticle, start_t_p, size_t_p, &zParticle[0]));
	NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zLevelParticle, start_t_p, size_t_p, &zLevelParticle[0]));

	int dimids_glCellIdx[] = {dimid_particles};
	size_t start_glCellIdx[1]={0}, count_glCellIdx[1]={nParticles};
	NC_SAFE_CALL( nc_get_vara_int(ncid, varid_glCellIdx, start_glCellIdx, count_glCellIdx, &glCellIdx[0]));

	NC_SAFE_CALL(nc_close(ncid));


	for (size_t i=0; i<glCellIdx.size(); i++){

		 // check if particle belongs to current block's gid
		 if (gcIdxToGid[glCellIdx[i]] == gid){
	 		EndPt p;
            p.pid = init;
            p.sid = init;
            p[0] = xParticle[i];  p[1] = yParticle[i];  p[2] = zParticle[i];
            p.zLevelParticle = zLevelParticle[i];
			p.glCellIdx = glCellIdx[i];
			particles.push_back(p);	
			init++;
		 }
	}

	dprint("particles %ld", particles.size());


}