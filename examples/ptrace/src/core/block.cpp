
#include "block.h"
#include "misc.h"
#include <netcdf.h>
#include "def.h"
#include <Eigen/Dense>
#include "nabo/nabo.h"
#include <map>
#include <utility>
#include <pnetcdf.h>

std::vector<Halo> block::get_halo_info()
{

	std::map<int, std::pair<std::vector<int>, std::vector<int>>> cells_verts;

	// iterate over cells
	for (size_t i = 0; i < indexToCellID.size(); i++)
	{

		int cellID = indexToCellID[i];
		// iterate over cell neighbors
		for (int j = i * maxEdges; j < (i + 1) * maxEdges; j++)
		{

			int nbr_cellID = cellsOnCell[j];
			// make sure nbr_cell is valid
			if (nbr_cellID != 0)
			{
				int nbr_gid = gcIdxToGid[nbr_cellID - 1];
				// if neighbor gid != current gid add then add cell id request to nbr_gid
				if (nbr_gid != gid)
				{
					cells_verts[nbr_gid].first.push_back(nbr_cellID);
				}
			}
		}

		// iterate over vertices on cell
		for (size_t j = i * maxEdges; j < (i + 1) * maxEdges; j++)
		{

			int nbr_vertexID = verticesOnCell[j];

			// make sure nbr_vertex is valid
			if (nbr_vertexID != 0)
			{
				int nbr_gid = gVIdxToGid[nbr_vertexID - 1];

				// if neighbor gid != current gid add then add vertex id request to nbr_gid
				if (nbr_gid != gid)
				{
					cells_verts[nbr_gid].second.push_back(nbr_vertexID);
				}
			}
		}
	}

	// iterate over map and populate halos

	std::vector<Halo> halos;

	std::for_each(cells_verts.begin(), cells_verts.end(),
				  [&](std::pair<int, std::pair<std::vector<int>, std::vector<int>>> element) {
					  Halo h;
					  h.src_gid = gid;
					  h.dest_gid = element.first;
					  h.glCellIDs = std::move(element.second.first);
					  h.glVertexIDs = std::move(element.second.second);

					  halos.push_back(h);
				  });

	return halos;
}

void block::process_halo_req(Halo &h, int framenum)
{

	// store the shell to populate dynamic halo containers later for velocityXv,...
	halo_info.push_back(h);

	// populate cell info for requested cellids
	for (size_t i = 0; i < h.glCellIDs.size(); i++)
	{
		h.xCell.push_back(xCell[cellIndex[h.glCellIDs[i]]]);
		h.yCell.push_back(yCell[cellIndex[h.glCellIDs[i]]]);
		h.zCell.push_back(zCell[cellIndex[h.glCellIDs[i]]]);
	}

	// populate vert info for requested vertex ids

	for (size_t i = 0; i < h.glVertexIDs.size(); i++)
	{

		h.xVertex.push_back(xVertex[vertexIndex[h.glVertexIDs[i]]]);
		h.yVertex.push_back(yVertex[vertexIndex[h.glVertexIDs[i]]]);
		h.zVertex.push_back(zVertex[vertexIndex[h.glVertexIDs[i]]]);
		for (size_t j = 0; j < nVertLevels; j++)
		{
			h.velocityXv.push_back(velocityXv[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]]);
			h.velocityYv.push_back(velocityYv[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]]);
			h.velocityZv.push_back(velocityZv[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]]);
			h.vertVelocityTop.push_back(vertVelocityTop[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]]);
			h.zMid.push_back(zMid[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]]);
			h.zTop.push_back(zTop[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]]);
		}
	}
}

void block::update_halo_info(Halo &h, int framenum)
{

	// iterate through each incoming cell and append info to local arrays

	size_t nCells_local = xCell.size();
	size_t nCells_all = nCells_local + h.glCellIDs.size();

	xCell.resize(nCells_all);
	yCell.resize(nCells_all);
	zCell.resize(nCells_all);

	for (size_t i = 0; i < h.glCellIDs.size(); i++)
	{

		int cellID = h.glCellIDs[i];

		cellIndex[cellID] = indexToCellID.size();
		indexToCellID.push_back(cellID);

		xCell[nCells_local + i] = h.xCell[i];
		yCell[nCells_local + i] = h.yCell[i];
		zCell[nCells_local + i] = h.zCell[i];
	}

	// iterate through each incoming vertex and append info to local arrays

	size_t nVertices_local = xVertex.size();
	size_t nVertices_all = nVertices_local + h.glVertexIDs.size();

	xVertex.resize(nVertices_all);
	yVertex.resize(nVertices_all);
	zVertex.resize(nVertices_all);
	velocityXv[framenum % 2].resize(nVertices_all * nVertLevels);
	velocityYv[framenum % 2].resize(nVertices_all * nVertLevels);
	velocityZv[framenum % 2].resize(nVertices_all * nVertLevels);
	vertVelocityTop[framenum % 2].resize(nVertices_all * nVertLevels);
	zMid[framenum % 2].resize(nVertices_all * nVertLevels);
	zTop[framenum % 2].resize(nVertices_all * nVertLevels);

	for (size_t i = 0; i < h.glVertexIDs.size(); i++)
	{

		int vertID = h.glVertexIDs[i];
		vertexIndex[vertID] = indexToVertexID.size();
		indexToVertexID.push_back(vertID);

		xVertex[nVertices_local + i] = h.xVertex[i];
		yVertex[nVertices_local + i] = h.yVertex[i];
		zVertex[nVertices_local + i] = h.zVertex[i];

		for (size_t j = 0; j < nVertLevels; j++)
		{
			velocityXv[framenum % 2][j + nVertLevels * (nVertices_local + i)] = h.velocityXv[j + nVertLevels * i];
			velocityYv[framenum % 2][j + nVertLevels * (nVertices_local + i)] = h.velocityYv[j + nVertLevels * i];
			velocityZv[framenum % 2][j + nVertLevels * (nVertices_local + i)] = h.velocityZv[j + nVertLevels * i];
			vertVelocityTop[framenum % 2][j + nVertLevels * (nVertices_local + i)] = h.vertVelocityTop[j + nVertLevels * i];
			zMid[framenum % 2][j + nVertLevels * (nVertices_local + i)] = h.zTop[j + nVertLevels * i];
			zTop[framenum % 2][j + nVertLevels * (nVertices_local + i)] = h.zTop[j + nVertLevels * i];
		}
	}
}

// prepare to send back dynamic request: velocityXv,... using halo_info
void block::process_halo_dynamic(int framenum)
{

	for (size_t g = 0; g < halo_info.size(); g++)
	{

		Halo &h = halo_info[g];
		size_t nVerts_to_send = h.glVertexIDs.size();

		h.velocityXv.resize(nVerts_to_send * nVertLevels);
		h.velocityYv.resize(nVerts_to_send * nVertLevels);
		h.velocityZv.resize(nVerts_to_send * nVertLevels);
		h.zMid.resize(nVerts_to_send * nVertLevels);
		h.zTop.resize(nVerts_to_send * nVertLevels);
		h.vertVelocityTop.resize(nVerts_to_send * nVertLevels);

		// iterate over vertices
		for (size_t i = 0; i < nVerts_to_send; i++)
		{

			for (size_t j = 0; j < nVertLevels; j++)
			{

				h.velocityXv[j] = velocityXv[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]];
				h.velocityYv[j] = velocityYv[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]];
				h.velocityZv[j] = velocityZv[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]];
				h.vertVelocityTop[j] = vertVelocityTop[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]];
				h.zMid[j] = zMid[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]];
				h.zTop[j] = zTop[framenum % 2][j + nVertLevels * vertexIndex[h.glVertexIDs[i]]];
			}
		}
	}
}

// do a dynamic update of velocityXv,...
void block::update_halo_dynamic(Halo &h, int framenum)
{

	size_t nVertices_all = indexToVertexID.size(); // should have been updated during static update rounds: during calls to update_halo_info()

	xVertex.resize(nVertices_all);
	yVertex.resize(nVertices_all);
	zVertex.resize(nVertices_all);
	velocityXv[framenum % 2].resize(nVertices_all * nVertLevels);
	velocityYv[framenum % 2].resize(nVertices_all * nVertLevels);
	velocityZv[framenum % 2].resize(nVertices_all * nVertLevels);
	vertVelocityTop[framenum % 2].resize(nVertices_all * nVertLevels);
	zMid[framenum % 2].resize(nVertices_all * nVertLevels);
	zTop[framenum % 2].resize(nVertices_all * nVertLevels);

	for (size_t i = 0; i < h.glVertexIDs.size(); i++)
	{

		int vertID = h.glVertexIDs[i];
		// vertexIndex[vertID] = indexToVertexID.size();
		// indexToVertexID.push_back(vertID);

		for (size_t j = 0; j < nVertLevels; j++)
		{

			velocityXv[framenum % 2][j + nVertLevels * vertexIndex[vertID]] = h.velocityXv[j + nVertLevels * i];
			velocityYv[framenum % 2][j + nVertLevels * vertexIndex[vertID]] = h.velocityYv[j + nVertLevels * i];
			velocityZv[framenum % 2][j + nVertLevels * vertexIndex[vertID]] = h.velocityZv[j + nVertLevels * i];
			zMid[framenum % 2][j + nVertLevels * vertexIndex[vertID]] = h.zMid[j + nVertLevels * i];
			zTop[framenum % 2][j + nVertLevels * vertexIndex[vertID]] = h.zTop[j + nVertLevels * i];
			vertVelocityTop[framenum % 2][j + nVertLevels * vertexIndex[vertID]] = h.vertVelocityTop[j + nVertLevels * i];
		}
	}
}

// creates the cellsOnCell from file (as opposed to from the incoming data)
void block::create_links(const std::string &fname_graph, const std::string &fname_graphpart, std::set<int> &links){

	std::vector<std::vector<int>> cell_nbrs = read_csv(fname_graph, ' ');
	// cellsOnCell.clear();
	// cellsOnCell.resize(maxEdges*indexToCellID.size());

	std::vector<std::vector<int>> partn_ids = read_csv(fname_graphpart.c_str());

	for (size_t i = 0; i < partn_ids.size(); i++)
		gcIdxToGid.push_back(partn_ids[i][0]);


	int idx = 0;
	for (size_t i = 0; i < partn_ids.size(); i++)
	{
		
		if (gid == partn_ids[i][0])
		{
			
			
			// iterate over neighbors
			for (int j=0; j<cell_nbrs[i+1].size(); j++){
				
				int nbr_cgid = cell_nbrs[i+1][j];
				// cellsOnCell[idx*maxEdges+j] = nbr_cgid;
				// check if neighbor in different partition
				if (partn_ids[nbr_cgid-1][0] != gid){
					links.insert(partn_ids[nbr_cgid-1][0]);
				}
			}
			idx++;

		}
	}


	

}

void block::create_links_mpas(const std::string &fname_graphinfo, std::set<int> &links, diy::mpi::communicator &world)
{

	std::vector<std::vector<int>> partn_ids = read_csv(fname_graphinfo.c_str());

	for (size_t i = 0; i < partn_ids.size(); i++)
		gcIdxToGid.push_back(partn_ids[i][0]);

	for (size_t i = 0; i < partn_ids.size(); i++)
	{
		if (gid == partn_ids[i][0])
		{
			for (int j = cellIndex[i + 1] * maxEdges; j < (cellIndex[i + 1] + 1) * maxEdges; j++)
			{
				if (cellsOnCell[j] != 0)
				{
					if (partn_ids[cellsOnCell[j] - 1][0] != gid)
					{
						links.insert(partn_ids[cellsOnCell[j] - 1][0]);
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
	for (size_t i = 0; i < all_indexToVertexID.size(); i++)
	{
		for (size_t j = 0; j < all_indexToVertexID[i].size(); j++)
		{

			gVIdxToGid[all_indexToVertexID[i][j] - 1] = all_gids[i];
		}
	}
}

// create a new particle file which includes a glCellIdx (global cell index) field
void block::generate_new_particle_file()
{

	if (gid == 0)
	{

		// read output.nc & create kd-tree

		int ncid;
		int dimid_cells;
		int varid_xCell, varid_yCell, varid_zCell;
		size_t nCells;

		std::vector<int> glCellIdx;

		// output2 to deal with online updation of output.nc by simulation
		NC_SAFE_CALL(nc_open("output2.nc", NC_NOWRITE, &ncid));

		NC_SAFE_CALL(nc_inq_dimid(ncid, "nCells", &dimid_cells));
		NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_cells, &nCells));

		NC_SAFE_CALL(nc_inq_varid(ncid, "xCell", &varid_xCell));
		NC_SAFE_CALL(nc_inq_varid(ncid, "yCell", &varid_yCell));
		NC_SAFE_CALL(nc_inq_varid(ncid, "zCell", &varid_zCell));

		dprint("nCells %ld", nCells);

		std::vector<double> xCells, yCells, zCells;
		xCells.resize(nCells);
		yCells.resize(nCells);
		zCells.resize(nCells);

		const size_t start_cells[1] = {0}, size_cells[1] = {nCells};
		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &xCells[0]));
		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &yCells[0]));
		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &zCells[0]));

		// dprint("varids %d %d %d,  %f %f %f", varid_xCell, varid_yCell, varid_zCell, xCells[1], yCells[1], zCells[1]);

		NC_SAFE_CALL(nc_close(ncid));

		/* KDTree */
		Eigen::MatrixXd C_local;
		Nabo::NNSearchD *nns_cells;

		C_local.resize(3, nCells);
		for (int i = 0; i < nCells; i++)
		{
			C_local(0, i) = xCells[i];
			C_local(1, i) = yCells[i];
			C_local(2, i) = zCells[i];
		}

		nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);

		// read particles.nc

		int dimid_particles, dimid_time;
		int varid_xParticle, varid_yParticle, varid_zParticle, varid_zLevelParticle;
		std::vector<double> xParticle, yParticle, zParticle, zLevelParticle;
		size_t nParticles;
		size_t nTime;

		NC_SAFE_CALL(nc_open("particles.nc", NC_CLOBBER, &ncid));

		NC_SAFE_CALL(nc_inq_dimid(ncid, "nParticles", &dimid_particles));
		NC_SAFE_CALL(nc_inq_dimid(ncid, "Time", &dimid_time));

		NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_particles, &nParticles));
		NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_time, &nTime));

		NC_SAFE_CALL(nc_inq_varid(ncid, "xParticle", &varid_xParticle));
		NC_SAFE_CALL(nc_inq_varid(ncid, "yParticle", &varid_yParticle));
		NC_SAFE_CALL(nc_inq_varid(ncid, "zParticle", &varid_zParticle));
		NC_SAFE_CALL(nc_inq_varid(ncid, "zLevelParticle", &varid_zLevelParticle));

		xParticle.resize(nTime * nParticles);
		yParticle.resize(nTime * nParticles);
		zParticle.resize(nTime * nParticles);
		zLevelParticle.resize(nTime * nParticles);
		glCellIdx.resize(nParticles);
		const size_t start_t_p[2] = {0, 0}, size_t_p[2] = {nTime, nParticles};

		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_xParticle, start_t_p, size_t_p, &xParticle[0]));
		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_yParticle, start_t_p, size_t_p, &yParticle[0]));
		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zParticle, start_t_p, size_t_p, &zParticle[0]));
		NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zLevelParticle, start_t_p, size_t_p, &zLevelParticle[0]));

		NC_SAFE_CALL(nc_close(ncid));

		// parce, filter and pushback particles

		int cellID;
		for (size_t i = 0; i < xParticle.size(); i++)
		{
			Eigen::VectorXd q(3);
			q << xParticle[i], yParticle[i], zParticle[i];
			// get nearest cell neighbor ids using mpas_c
			Eigen::VectorXi nearest_cell_idx(1);
			Eigen::VectorXd dists2_cell(1);
			nns_cells->knn(q, nearest_cell_idx, dists2_cell, 1, 0, Nabo::NNSearchF::SORT_RESULTS | Nabo::NNSearchF::ALLOW_SELF_MATCH);
			cellID = nearest_cell_idx[0]; // cell global idx
			glCellIdx[i] = cellID + 1;
			dprint("i %ld cellID %d, %f %f %f, %f %f %f", i, cellID, xCells[i], yCells[i], zCells[i], xParticle[i], yParticle[i], zParticle[i]);
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
		start[0] = 0;
		start[1] = 0;
		count[0] = nTime;
		count[1] = nParticles;
		NC_SAFE_CALL(nc_put_vara_double(ncid, varid_xParticle, start, count, &xParticle[0]));
		NC_SAFE_CALL(nc_put_vara_double(ncid, varid_yParticle, start, count, &yParticle[0]));
		NC_SAFE_CALL(nc_put_vara_double(ncid, varid_zParticle, start, count, &zParticle[0]));
		NC_SAFE_CALL(nc_put_vara_double(ncid, varid_zLevelParticle, start, count, &zLevelParticle[0]));

		size_t start_glCellIdx[1] = {0}, count_glCellIdx[1] = {nParticles};
		NC_SAFE_CALL(nc_put_vara_int(ncid, varid_glCellIdx, start_glCellIdx, count_glCellIdx, &glCellIdx[0]));

		NC_SAFE_CALL(nc_close(ncid));

	} // if (gid==0)
}

void block::init_seeds_particles(diy::mpi::communicator& world, std::string &fname_particles, int framenum){

	// initialize seeds directly from particles.nc

	int ncid;
	int dimid_particles, dimid_time;
	int varid_xParticle, varid_yParticle, varid_zParticle, varid_zLevelParticle, varid_glCellIdx;
	std::vector<double> xParticle, yParticle, zParticle, zLevelParticle;
	std::vector<int> glCellIdx;
	std::vector<int> currentBlock;
	int varid_currentBlock;

	MPI_Offset nParticles;
	MPI_Offset nTime;

	PNC_SAFE_CALL(ncmpi_open(world, fname_particles.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid));

	PNC_SAFE_CALL(ncmpi_inq_dimid(ncid, "nParticles", &dimid_particles));
	PNC_SAFE_CALL(ncmpi_inq_dimid(ncid, "Time", &dimid_time));

	PNC_SAFE_CALL(ncmpi_inq_dimlen(ncid, dimid_particles, &nParticles));
	PNC_SAFE_CALL(ncmpi_inq_dimlen(ncid, dimid_time, &nTime));

	PNC_SAFE_CALL(ncmpi_inq_varid(ncid, "xParticle", &varid_xParticle));
	PNC_SAFE_CALL(ncmpi_inq_varid(ncid, "yParticle", &varid_yParticle));
	PNC_SAFE_CALL(ncmpi_inq_varid(ncid, "zParticle", &varid_zParticle));
	PNC_SAFE_CALL(ncmpi_inq_varid(ncid, "zLevelParticle", &varid_zLevelParticle));
	PNC_SAFE_CALL(ncmpi_inq_varid(ncid, "currentCell", &varid_glCellIdx));

	xParticle.resize(nTime * nParticles);
	yParticle.resize(nTime * nParticles);
	zParticle.resize(nTime * nParticles);
	zLevelParticle.resize(nTime * nParticles);
	glCellIdx.resize(nTime * nParticles);
	currentBlock.resize(nTime * nParticles);
	const MPI_Offset start_t_p[2] = {0, 0}, size_t_p[2] = {nTime, nParticles};

	PNC_SAFE_CALL(ncmpi_get_vara_double_all(ncid, varid_xParticle, start_t_p, size_t_p, &xParticle[0]));
	PNC_SAFE_CALL(ncmpi_get_vara_double_all(ncid, varid_yParticle, start_t_p, size_t_p, &yParticle[0]));
	PNC_SAFE_CALL(ncmpi_get_vara_double_all(ncid, varid_zParticle, start_t_p, size_t_p, &zParticle[0]));
	PNC_SAFE_CALL(ncmpi_get_vara_double_all(ncid, varid_zLevelParticle, start_t_p, size_t_p, &zLevelParticle[0]));

	PNC_SAFE_CALL(ncmpi_get_vara_int_all(ncid, varid_glCellIdx, start_t_p, size_t_p, &glCellIdx[0]));


	// int ncid_p;

	// NC_SAFE_CALL(nc_open("particles.nc", NC_CLOBBER, &ncid_p));
	PNC_SAFE_CALL(ncmpi_inq_varid(ncid, "currentBlock", &varid_currentBlock));
	PNC_SAFE_CALL(ncmpi_get_vara_int_all(ncid, varid_currentBlock, start_t_p, size_t_p, &currentBlock[0]));
	PNC_SAFE_CALL(ncmpi_close(ncid));	



	for (size_t i=0; i<xParticle.size(); i++){
		if (world.rank() == currentBlock[i] && init<2000 ){//&& init == 37){
			EndPt p;
			p.pid = init;
			p.sid = init;
			p[0] = xParticle[i];
			p[1] = yParticle[i];
			p[2] = zParticle[i];
			if (pow(-6452390 - p[0], 2)+ pow(-447879 - p[1], 2) + pow(197714 - p[2],2)<pow(5524190.98744,2)){
				p.zLevelParticle = zLevelParticle[i];
				p.glCellIdx = glCellIdx[i];
				particles.push_back(p);
				dprint("Init %d cellid %d", init, p.glCellIdx);
			}
		
		}
		init++;	

	}
	dprint("particle %ld", particles.size());


}

void block::init_seeds_mpas(std::string &fname_particles, int framenum, int rank)
{

	// read particles2.nc
	int ncid;
	int dimid_particles, dimid_time;
	int varid_xParticle, varid_yParticle, varid_zParticle, varid_zLevelParticle, varid_glCellIdx;
	std::vector<double> xParticle, yParticle, zParticle, zLevelParticle;
	std::vector<int> glCellIdx;
	std::vector<int> currentBlock;
	int varid_currentBlock;

	size_t nParticles;
	size_t nTime;

	NC_SAFE_CALL(nc_open(fname_particles.c_str(), NC_CLOBBER, &ncid));

	NC_SAFE_CALL(nc_inq_dimid(ncid, "nParticles", &dimid_particles));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "Time", &dimid_time));

	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_particles, &nParticles));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_time, &nTime));

	NC_SAFE_CALL(nc_inq_varid(ncid, "xParticle", &varid_xParticle));
	NC_SAFE_CALL(nc_inq_varid(ncid, "yParticle", &varid_yParticle));
	NC_SAFE_CALL(nc_inq_varid(ncid, "zParticle", &varid_zParticle));
	NC_SAFE_CALL(nc_inq_varid(ncid, "zLevelParticle", &varid_zLevelParticle));
	NC_SAFE_CALL(nc_inq_varid(ncid, "glCellIdx", &varid_glCellIdx));

	xParticle.resize(nTime * nParticles);
	yParticle.resize(nTime * nParticles);
	zParticle.resize(nTime * nParticles);
	zLevelParticle.resize(nTime * nParticles);
	glCellIdx.resize(nParticles);
	currentBlock.resize(nTime * nParticles);
	const size_t start_t_p[2] = {0, 0}, size_t_p[2] = {nTime, nParticles};

	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_xParticle, start_t_p, size_t_p, &xParticle[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_yParticle, start_t_p, size_t_p, &yParticle[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zParticle, start_t_p, size_t_p, &zParticle[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zLevelParticle, start_t_p, size_t_p, &zLevelParticle[0]));

	int dimids_glCellIdx[] = {dimid_particles};
	size_t start_glCellIdx[1] = {0}, count_glCellIdx[1] = {nParticles};
	NC_SAFE_CALL(nc_get_vara_int(ncid, varid_glCellIdx, start_glCellIdx, count_glCellIdx, &glCellIdx[0]));

	NC_SAFE_CALL(nc_close(ncid));

	int ncid_p;

	NC_SAFE_CALL(nc_open("particles.nc", NC_CLOBBER, &ncid_p));
	NC_SAFE_CALL(nc_inq_varid(ncid, "currentBlock", &varid_currentBlock));
	NC_SAFE_CALL(nc_get_vara_int(ncid, varid_currentBlock, start_t_p, size_t_p, &currentBlock[0]));
	NC_SAFE_CALL(nc_close(ncid_p));

	for (size_t i = 0; i < glCellIdx.size(); i++)
	{

		// check if particle belongs to current block's gid
		//  if (gcIdxToGid[glCellIdx[i]] == gid){

		// if (i%32==0) continue;

		// if (init==130){
		// if (init==410 || init == 130){
        
		if (rank == currentBlock[i]){
			EndPt p;
			p.pid = init;
			p.sid = init;
			p[0] = xParticle[i];
			p[1] = yParticle[i];
			p[2] = zParticle[i];
			p.zLevelParticle = zLevelParticle[i];
			p.glCellIdx = glCellIdx[i];
			particles.push_back(p);
		// dprint("Init cellid %d", p.glCellIdx);
		}
		init++;

		// if (init>230) break;
		//  }
	}

	// particles.resize(1);
	dprint("particles %ld, framenum %d, rank %d", particles.size(), framenum, rank);
}

static void handle_error(int status, int lineno)
{
	fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
	MPI_Abort(MPI_COMM_WORLD, 1);
}

// void block::parallel_write_segments(MPI_Comm comm, int max_steps){
void block::parallel_write_segments(diy::mpi::communicator &comm, int max_steps)
{

	int ret, ncfile, nprocs, rank;
	char filename[256];
	int data[2];

	// MPI_Comm_rank(comm, &rank);
	// MPI_Comm_size(comm, &nprocs);

	rank = comm.rank();
	nprocs = comm.size();

	strcpy(filename, "segments.nc");

	std::vector<int> segsizes(nprocs); // later contain total length of all segments in each proc
	std::vector<int> segsizes_local;   // sizes of local segments
	// std::vector<int> segstep_local; // step values of local segments
	std::vector<int> segsid_local; // seeds of local segments
	int totalseglen_local = 0;	 // later contain total length of all local segments

	for (size_t i = 0; i < segments.size(); i++)
	{
		totalseglen_local += segments[i].pts.size();
		segsizes_local.push_back(segments[i].pts.size());
		// segstep_local.push_back(segments[i].step);
		// segsid_local.push_back(segments[i].sid);
		segsid_local.push_back(segments[i].pid);
	}

	// for (size_t i=0; i<segsizes_local.size(); i++){
	// 	dprint("rank %d i %ld segl %d", rank, i, segsizes_local[i]);
	// }

	// dprint("rank %d, totalseglen_local %d", rank, totalseglen_local);

	long nSegments_global;
	long nSegments_local = segments.size();
	// get num of segments globally
	MPI_Allreduce(
		&nSegments_local,
		&nSegments_global,
		1,
		MPI_LONG,
		MPI_SUM,
		comm);

	long nSegments_global_granular[nprocs];
	MPI_Allgather(
		&nSegments_local,
		1,
		MPI_LONG,
		&nSegments_global_granular[0],
		1,
		MPI_LONG,
		comm);

	std::vector<int> segcnt_offsets(nprocs); // segment count offsets
	segcnt_offsets[0] = 0;
	long max_nSegments = nSegments_global_granular[0]; // later has maximum number of segments in any process
	for (size_t i = 1; i < nprocs; i++)
	{
		segcnt_offsets[i] = segcnt_offsets[i - 1] + nSegments_global_granular[i - 1];
		if (max_nSegments < nSegments_global_granular[i])
			max_nSegments = nSegments_global_granular[i];
	}
	// for (int i=0; i<nprocs;i++){
	// 	dprint("i %d nsg %d", i, segcnt_offsets[i]);
	// }

	// get number of steps globally split by process into segsizes
	MPI_Allgather(
		&totalseglen_local,
		1,
		MPI_INT,
		&segsizes[0],
		1,
		MPI_INT,
		comm);

	std::vector<int> seg_offsets(nprocs);
	seg_offsets[0] = 0;
	int totalseglen_global = segsizes[0];
	for (size_t i = 1; i < segsizes.size(); i++)
	{
		seg_offsets[i] = seg_offsets[i - 1] + segsizes[i - 1];
		totalseglen_global += segsizes[i];
	}
	dprint("totalseglen_global %d", totalseglen_global);
	// dprint("nSegments_global %d", nSegments_global);

	// if (rank==7)
	// for (size_t i=0; i<segments.size(); i++){
	// 	dprint("rank %d i %ld seg %ld", rank, i, segments[i].pts.size());
	// }

	// dprint("rank %d segsizes size %d %d %d %d %d %d %d %d %ld", rank, segsizes[0], segsizes[1], segsizes[2], segsizes[3], segsizes[4], segsizes[5], segsizes[6], segsizes[7], segsizes.size());

	ret = ncmpi_create(comm, filename,
					   NC_CLOBBER | NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	int dimid_x, dimid_y, dimid_z, dimid_segsize;
	int varid_segsizes, varid_seedid, varid_step, varid_px, varid_py, varid_pz;

	// define dimensions

	ret = ncmpi_def_dim(ncfile, "x", totalseglen_global, &dimid_x);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	ret = ncmpi_def_dim(ncfile, "y", totalseglen_global, &dimid_y);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	ret = ncmpi_def_dim(ncfile, "z", totalseglen_global, &dimid_z);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	ret = ncmpi_def_dim(ncfile, "segsize", nSegments_global, &dimid_segsize);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	// define variables
	int ndims_pos = 1;
	int pos_dims[1] = {dimid_x};
	ret = ncmpi_def_var(ncfile, "px", NC_DOUBLE, ndims_pos, pos_dims, &varid_px);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	pos_dims[0] = {dimid_y};
	ret = ncmpi_def_var(ncfile, "py", NC_DOUBLE, ndims_pos, pos_dims, &varid_py);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	pos_dims[0] = {dimid_z};
	ret = ncmpi_def_var(ncfile, "pz", NC_DOUBLE, ndims_pos, pos_dims, &varid_pz);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	pos_dims[0] = {dimid_segsize};
	ret = ncmpi_def_var(ncfile, "segsizes", NC_INT, ndims_pos, pos_dims, &varid_segsizes);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	pos_dims[0] = {dimid_segsize};
	ret = ncmpi_def_var(ncfile, "seedid", NC_INT, ndims_pos, pos_dims, &varid_seedid);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	pos_dims[0] = {dimid_segsize};
	ret = ncmpi_def_var(ncfile, "step", NC_INT, ndims_pos, pos_dims, &varid_step);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	// end define mode
	ret = ncmpi_enddef(ncfile);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	// allocate starts and counts, allocate buffer, write
	// dprint("rank %d segoffset[rank] %d", rank, seg_offsets[rank]);
	MPI_Offset start_segsize[1] = {static_cast<long long int>(segcnt_offsets[rank])}, count_segsize[1] = {static_cast<long long int>(segsizes_local.size())};

	dprint("rank %d segsizes start %lld count %ld, nseg_glob %ld", comm.rank(), segcnt_offsets[rank], segsizes_local.size(), nSegments_global);
	// if (segcnt_offsets[rank] == nSegments_global && segsizes_local.size() == 0)
	// 	ret = ncmpi_put_vara_int_all(ncfile, varid_segsizes, start_segsize, count_segsize, &segsizes_local[0]);
	// else 
		ret = ncmpi_put_vara_int_all(ncfile, varid_segsizes, start_segsize, count_segsize, &segsizes_local[0]);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);
	dprint("first int");
	// ret = ncmpi_put_vara_int_all(ncfile, varid_step, start_segsize, count_segsize, &segstep_local[0]);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);
	dprint("second int");
	ret = ncmpi_put_vara_int_all(ncfile, varid_seedid, start_segsize, count_segsize, &segsid_local[0]);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);

	dprint("here rank %d segsize %ld", rank, segments.size());
	long long int local_offset = 0;
	for (size_t i = 0; i < segments.size(); i++)
	{

		MPI_Offset start[1] = {seg_offsets[rank] + local_offset}, count[1] = {(long long int)segments[i].pts.size()};

		double buffer_x[segments[i].pts.size()];
		double buffer_y[segments[i].pts.size()];
		double buffer_z[segments[i].pts.size()];

		for (size_t j = 0; j < segments[i].pts.size(); j++)
		{
			buffer_x[j] = segments[i].pts[j].coords[0];
			buffer_y[j] = segments[i].pts[j].coords[1];
			buffer_z[j] = segments[i].pts[j].coords[2];
		}
		// if (rank==5)
		// 	dprint("rank %d, segsize %ld, i %ld, start %lld, size %ld", rank, segments[i].pts.size(), i, seg_offsets[rank]+local_offset, segments[i].pts.size());

		ret = ncmpi_put_vara_double_all(ncfile, varid_px, start, count, &buffer_x[0]);
		if (ret != NC_NOERR)
			handle_error(ret, __LINE__);

		ret = ncmpi_put_vara_double_all(ncfile, varid_py, start, count, &buffer_y[0]);
		if (ret != NC_NOERR)
			handle_error(ret, __LINE__);

		ret = ncmpi_put_vara_double_all(ncfile, varid_pz, start, count, &buffer_z[0]);
		if (ret != NC_NOERR)
			handle_error(ret, __LINE__);

		local_offset += segments[i].pts.size();
	}

	size_t extra_runs = max_nSegments - segments.size();
	for (size_t i = 0; i < extra_runs; i++)
	{

		MPI_Offset start[1] = {0}, count[1] = {0};
		double buffer_x[1];
		double buffer_y[1];
		double buffer_z[1];

		ret = ncmpi_put_vara_double_all(ncfile, varid_px, start, count, &buffer_x[0]);
		if (ret != NC_NOERR)
			handle_error(ret, __LINE__);

		ret = ncmpi_put_vara_double_all(ncfile, varid_py, start, count, &buffer_y[0]);
		if (ret != NC_NOERR)
			handle_error(ret, __LINE__);

		ret = ncmpi_put_vara_double_all(ncfile, varid_pz, start, count, &buffer_z[0]);
		if (ret != NC_NOERR)
			handle_error(ret, __LINE__);
	}

	// close file
	ret = ncmpi_close(ncfile);
	if (ret != NC_NOERR)
		handle_error(ret, __LINE__);
}

void block::parallel_write_simstep_segments(diy::mpi::communicator &comm, int framenum)
{

	// // 	dprint("new writer framenum %d", framenum);

	// int ret, ncfile, nprocs, rank;
	// char filename[256];

	// rank = comm.rank();
	// nprocs = comm.size();

	// strcpy(filename, "segments2.nc");

	// // create if frame==2 create; afterward open and update file

	// int dimid_nParticles, dimid_Time, dimid_Epochs;

	// int varid_x, varid_y, varid_z, varid_indexToPid;
	// int ndims = 2;
	// // dprint("nParticles %ld", particles.size());
	// if (framenum == 2)
	// {
	// 	size_t nParticles, Time;
	// 	NC_SAFE_CALL(ncmpi_create(comm, filename, NC_CLOBBER, MPI_INFO_NULL, &ncfile));

	// 	diy::mpi::all_reduce(comm, particles.size(), nParticles, std::plus<size_t>());

	// 	// nParticles = 164;

	// 	NC_SAFE_CALL(ncmpi_def_dim(ncfile, "nParticles", nParticles, &dimid_nParticles));
	// 	NC_SAFE_CALL(ncmpi_def_dim(ncfile, "Time", NC_UNLIMITED, &dimid_Time));

	// 	int pos_dims[] = {dimid_Time, dimid_nParticles};
	// 	NC_SAFE_CALL(ncmpi_def_var(ncfile, "x", NC_DOUBLE, ndims, pos_dims, &varid_x));
	// 	NC_SAFE_CALL(ncmpi_def_var(ncfile, "y", NC_DOUBLE, ndims, pos_dims, &varid_y));
	// 	NC_SAFE_CALL(ncmpi_def_var(ncfile, "z", NC_DOUBLE, ndims, pos_dims, &varid_z));

	// 	// int indexToPid_dims [] = {dimid_nParticles};
	// 	NC_SAFE_CALL(ncmpi_def_var(ncfile, "indexToPid", NC_INT, ndims, pos_dims, &varid_indexToPid));

	// 	// end define mode
	// 	NC_SAFE_CALL(ncmpi_enddef(ncfile));

	// 	// close file
	// 	NC_SAFE_CALL(ncmpi_close(ncfile));
	// }

	// size_t start_idx = 0;
	// diy::mpi::scan(comm, segments.size(), start_idx, std::plus<size_t>());
	// start_idx -= segments.size();
	// // dprint("here %ld", start_idx);

	// int nSegments_local = segments.size();

	// ncfile = 0;
	// NC_SAFE_CALL(ncmpi_open(comm, filename, NC_WRITE, MPI_INFO_NULL, &ncfile));
	// int cur_segsize = 0;
	// if (nSegments_local > 0)
	// {
	// 	cur_segsize = (int)segments[0].pts.size(); // number of time steps per segment
	// }

	// long long int nParticles;
	// NC_SAFE_CALL(ncmpi_inq_dimid(ncfile, "nParticles", &dimid_nParticles));
	// NC_SAFE_CALL(ncmpi_inq_dimlen(ncfile, dimid_nParticles, &nParticles));

	// long long int Time;
	// NC_SAFE_CALL(ncmpi_inq_dimid(ncfile, "Time", &dimid_Time));
	// NC_SAFE_CALL(ncmpi_inq_dimlen(ncfile, dimid_Time, &Time));
	// NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "x", &varid_x));
	// NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "y", &varid_y));
	// NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "z", &varid_z));
	// NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "indexToPid", &varid_indexToPid));

	// MPI_Offset start[ndims], count[ndims];
	// start[0] = (framenum - 2) * cur_segsize;
	// count[0] = cur_segsize;
	// start[1] = start_idx;
	// count[1] = segments.size();

	// //  write segment data for current round and clear segments

	// double *data_x = NULL, *data_y = NULL, *data_z = NULL;
	// std::vector<int> indexToPid(nSegments_local);
	// // 	int *cur_pids = NULL;

	// if (nSegments_local > 0)
	// {
	// 	data_x = new double[nSegments_local * cur_segsize];
	// 	data_y = new double[nSegments_local * cur_segsize];
	// 	data_z = new double[nSegments_local * cur_segsize];
	// 	// cur_pids = new int[cur_nSegments];
	// 	for (size_t i = 0; i < segments.size(); i++)
	// 	{
	// 		for (size_t j = 0; j < segments[i].pts.size(); j++)
	// 		{

	// 			data_x[i + j * nSegments_local] = segments[i].pts[j].coords[0];
	// 			data_y[i + j * nSegments_local] = segments[i].pts[j].coords[1];
	// 			data_z[i + j * nSegments_local] = segments[i].pts[j].coords[2];
	// 		}
	// 		// cur_pids[i] = segments[i].pid;
	// 	}

	// 	for (size_t i = 0; i < segments.size(); i++)
	// 	{
	// 		indexToPid[i] = segments[i].pid;
	// 	}
	// }

	// NC_SAFE_CALL(ncmpi_put_vara_double_all(ncfile, varid_x, &start[0], &count[0], data_x));
	// NC_SAFE_CALL(ncmpi_put_vara_double_all(ncfile, varid_y, &start[0], &count[0], data_y));
	// NC_SAFE_CALL(ncmpi_put_vara_double_all(ncfile, varid_z, &start[0], &count[0], data_z));

	// start[0] = (framenum - 2);
	// count[0] = 1;
	// start[1] = start_idx;
	// count[1] = segments.size();

	// NC_SAFE_CALL(ncmpi_put_vara_int_all(ncfile,varid_indexToPid, &start[0], &count[0], &indexToPid[0]));


	// // close file
	// NC_SAFE_CALL(ncmpi_close(ncfile));

	// segments.clear();

	// if (data_x != NULL)
	// 	delete[] data_x;
	// if (data_y != NULL)
	// 	delete[] data_y;
	// if (data_z != NULL)
	// 	delete[] data_z;
}

// void block::parallel_write_simstep_segments(diy::mpi::communicator &comm, int framenum){

// 	dprint("new writer framenum %d", framenum);

// 	int ret, ncfile, nprocs, rank;
// 	char filename[256];

// 	rank = comm.rank();
// 	nprocs = comm.size();

// 	strcpy(filename, "segments2.nc");

// 	// create if frame==2 create; afterward open and update file

// 	int dimid_nParticles, dimid_Time;

// 	int varid_x, varid_y, varid_z;
// 	int ndims = 2;
// 	// dprint("nParticles %ld", particles.size());
// 	if (framenum==2){
// 		size_t nParticles, Time;
// 		NC_SAFE_CALL(ncmpi_create(comm, filename, NC_CLOBBER, MPI_INFO_NULL, &ncfile));

// 		// diy::mpi::all_reduce(comm, particles.size(), nParticles, std::plus<size_t>());

// 		nParticles = 164;

// 		NC_SAFE_CALL(ncmpi_def_dim(ncfile, "nParticles", nParticles, &dimid_nParticles));
// 		NC_SAFE_CALL(ncmpi_def_dim(ncfile, "Time", NC_UNLIMITED, &dimid_Time));

// 		int pos_dims [] = {dimid_Time, dimid_nParticles};
// 		NC_SAFE_CALL(ncmpi_def_var(ncfile, "x", NC_DOUBLE, ndims, pos_dims, &varid_x));
// 		NC_SAFE_CALL(ncmpi_def_var(ncfile, "y", NC_DOUBLE, ndims, pos_dims, &varid_y));
// 		NC_SAFE_CALL(ncmpi_def_var(ncfile, "z", NC_DOUBLE, ndims, pos_dims, &varid_z));

// 		// end define mode
// 		NC_SAFE_CALL(ncmpi_enddef(ncfile));

// 		// close file
// 		NC_SAFE_CALL(ncmpi_close(ncfile));
// 	}

// 	dprint("segments size %ld", segments.size());

// 	ncfile = 0;
// 	NC_SAFE_CALL( ncmpi_open(comm, filename, NC_WRITE, MPI_INFO_NULL, &ncfile) );

// 	long long int nParticles;
// 	NC_SAFE_CALL(ncmpi_inq_dimid(ncfile, "nParticles", &dimid_nParticles) );
// 	NC_SAFE_CALL(ncmpi_inq_dimlen(ncfile, dimid_nParticles, &nParticles));

// 	long long int Time;
// 	NC_SAFE_CALL(ncmpi_inq_dimid(ncfile, "Time", &dimid_Time) );
// 	NC_SAFE_CALL(ncmpi_inq_dimlen(ncfile, dimid_Time, &Time));
// 	NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "x", &varid_x));
// 	NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "y", &varid_y));
// 	NC_SAFE_CALL(ncmpi_inq_varid(ncfile, "z", &varid_z));

// 	// write segment data for current round and clear segments

// 	dprint("one");

// 	size_t cur_nSegments = segments.size();
// 	double *data_x = NULL, *data_y = NULL, *data_z = NULL;
// 	size_t cur_segsize = 0;
// 	int *cur_pids = NULL;

// 	if (cur_nSegments > 0){
// 		cur_segsize = segments[0].pts.size(); // number of time steps per segment
// 		// data_x = new double[cur_nSegments*cur_segsize];
// 		// data_y = new double[cur_nSegments*cur_segsize];
// 		// data_z = new double[cur_nSegments*cur_segsize];
// 		cur_pids = new int[cur_nSegments];
// 		for (size_t i=0; i<segments.size(); i++){
// 			for (size_t j=0; j<segments[i].pts.size(); j++){
// 				// data_x[i*cur_segsize + j] = segments[i].pts[j].coords[0];
// 				// data_y[i*cur_segsize + j] = segments[i].pts[j].coords[1];
// 				// data_z[i*cur_segsize + j] = segments[i].pts[j].coords[2];

// 				// data_x[i + j*cur_nSegments] = segments[i].pts[j].coords[0];
// 				// data_y[i*cur_segsize + j] = segments[i].pts[j].coords[1];
// 				// data_z[i*cur_segsize + j] = segments[i].pts[j].coords[2];
// 			}
// 			cur_pids[i] = segments[i].pid;
// 		}

// 	}
// 	dprint("two, cur_nSegments %ld", cur_nSegments);

// 	double data[] = {0.1, 0.2, 0.3, 0.4};
// 	int num_reqs = (int) cur_nSegments; // set to number of points in current epoch in the current process, should be cur_nSegments

// 	// MPI_Offset start[num_reqs][ndims], count[num_reqs][ndims];
// 	MPI_Offset **starts = NULL, **counts = NULL;

// 	if (num_reqs > 0) {
// 		starts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
// 		starts[0] = (MPI_Offset*)  calloc(num_reqs*ndims, sizeof(MPI_Offset));
// 		for (int i=1; i<num_reqs; i++)
// 			starts[i] = starts[i-1] + ndims;

// 		counts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
// 		counts[0] = (MPI_Offset*)  calloc(num_reqs*ndims, sizeof(MPI_Offset));
// 		for (int i=1; i<num_reqs; i++)
// 			counts[i] = counts[i-1] + ndims;

// 	}

// 	dprint("three num_reqs %d, ndims %d", num_reqs, ndims);
// 	// counts[0][0] = 4;
// 	// counts[0][1] = 1;
// 	// starts[0][0] = (framenum-2) * 4;
// 	// starts[0][1] = framenum%2;

// 	if (cur_nSegments > 0)
// 	{
// 		for (int i = 0; i < num_reqs; i++)
// 		{
// 			// for (int j=0; j<ndims; j++){
// 			starts[i][0] = (framenum - 2) * cur_segsize;
// 			counts[i][0] = cur_segsize;
// 			starts[i][1] = cur_pids[i];
// 			counts[i][1] = 1;
// 			// }
// 			// dprint("%lld %lld %lld %lld", starts[i][0], starts[i][1], counts[i][0], counts[i][1]);
// 		}

// 		int buf_len = 0;
// 		for (int i = 0; i < num_reqs; i++)
// 		{
// 			MPI_Offset w_req_len = 1;
// 			for (int j = 0; j < ndims; j++)
// 				w_req_len *= counts[i][j];
// 			buf_len += w_req_len;
// 		}
// 		data_x = (double*) malloc(buf_len * sizeof(double));
// 		for (int i = 0; i < buf_len; i++)
// 			data_x[i] = (double)rank;

// 		dprint("four num_reqs %d", num_reqs);
// 		NC_SAFE_CALL(ncmpi_put_varn_double_all(ncfile, varid_x, num_reqs, starts, counts, data_x));
// 		// NC_SAFE_CALL(ncmpi_put_varn_double_all(ncfile, varid_y, num_reqs, starts, counts, &data[0]));
// 		// NC_SAFE_CALL(ncmpi_put_varn_double_all(ncfile, varid_z, num_reqs, starts, counts, &data[0]));
// 	}
// 	dprint("five");

// 	if (data_x != NULL) delete [] data_x; data_x = NULL;
// 	if (data_y != NULL) delete [] data_y; data_y = NULL;
// 	if (data_z != NULL) delete [] data_z; data_z = NULL;
// 	if (cur_pids != NULL) delete [] cur_pids; cur_pids = NULL;

// 	if (starts != NULL ) {
// 		if (starts[0] != NULL) delete [] starts[0]; // this if condition is moot
// 		starts[0] = NULL;
// 		delete [] starts;
// 		starts = NULL;

// 	}
// 	if (counts != NULL ){
// 		if (counts[0] != NULL) delete [] counts[0]; // this if contition is moot
// 		counts[0] = NULL;
// 		delete [] counts;
// 		counts = NULL;
// 	}
// 	dprint("six");
// 	// close file
// 	NC_SAFE_CALL(ncmpi_close(ncfile));

// 	segments.clear();

// }

//http://cucis.ece.northwestern.edu/projects/PnetCDF/doc/pnetcdf-c/ncmpi_005fput_005fvarn_005f_003ctype_003e.html#ncmpi_005fput_005fvarn_005f_003ctype_003e