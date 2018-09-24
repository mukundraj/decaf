#include <decaf/decaf.hpp>
#include <bredala/data_model/simplefield.hpp>
#include <bredala/data_model/vectorfield.hpp>
#include <bredala/data_model/boost_macros.h>

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi/datatypes.hpp>
#include <diy/io/bov.hpp>
#include <diy/pick.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <map>
#include <cstdlib>

using namespace decaf;
using namespace std;



#include<iostream>
#include<string>
#include "streamlines.h"
#include "pathlines.h"
#include "mpaso.h"
#include "pblock.hpp"
#include "block_io.h"

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/serialization.hpp>

#include <Eigen/Dense>
#include "nabo/nabo.h"
#include "interpolators.h"

//#include "advect.h"


  void parallel_streamlines(mpaso &mpas_g, mpaso &mpas_c, int round, PBlock &b, int max_steps, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner, diy::mpi::communicator &world );
int frame_number = 0, global_gid; //FIX_ME: remove global variable

void populate_block( PBlock* b,
		int gid,
		std::vector<double> &data_bar, mpaso *mpasobj, diy::mpi::communicator &world, mpaso &mpas_c ){


	int frame_no = (int) data_bar[13];
	int nVertices = (int) data_bar[14];
	int nCells = (int) data_bar[15];
	int nVertLevels = (int) data_bar[16];
	// int nCells_global;
	std::vector<int> nCell_counts, nVertices_counts; 
	std::vector<std::vector<double>> xCells_global, yCells_global, zCells_global, 
		xVertex_global, yVertex_global, zVertex_global; 

	int read_order[13] = {4, 5, 6, 11, 0, 1, 2, 3, 8, 9,10,7, 12};
	int field_offsets[13];
	int field_sizes[13];
	int bsp = 0; // block start position

	field_offsets[read_order[0]] = bsp+20;
	field_sizes[read_order[0]] = data_bar[bsp+read_order[0]];
	for (int j=1; j<13; j++){
		field_offsets[read_order[j]] = field_offsets[read_order[j-1]]+field_sizes[read_order[j-1]];
		field_sizes[read_order[j]] = data_bar[bsp+read_order[j]];
	}

	fprintf(stderr, "inside populate %d %d %d %d %d\n", gid, frame_no, nVertices, nCells, nVertices);
	if (frame_no==1){
		b->indexToCellID.resize(nCells);  //0  
		b->xCell.resize(nCells);  //1
		b->yCell.resize(nCells); //2
		b->zCell.resize(nCells); //3

		b->velocityX.resize(nVertLevels*nCells); //4
		b->velocityY.resize(nVertLevels*nCells); //5
		b->velocityZ.resize(nVertLevels*nCells); //6

		b->indexToVertexID.resize(nVertices);	//7
		b->xVertex.resize(nVertices); //8
		b->yVertex.resize(nVertices); //9
		b->zVertex.resize(nVertices); //10
		b->zTop.resize(nVertLevels*nCells); //11
		b->cellsOnVertex.resize(nVertices*3); //12

		// for (int i =0;i<13; i++)
		//   fprintf(stderr, "field offsets %d %d %d\n", i, field_offsets[i], field_sizes[i]);
		// populate frame 1 exclusive members 
		int k = 0; // xCell index offset var
		for (int j=field_offsets[0]; j<field_offsets[0]+field_sizes[0]; j++){
			b->xCell[k] = data_bar[field_offsets[1]+k]; // -1 on LHS to deal with fortran to C offset
			b->yCell[k] = data_bar[field_offsets[2]+k];
			b->zCell[k] = data_bar[field_offsets[3]+k];
			b->indexToCellID[k] = (int)data_bar[field_offsets[0]+k];
			// b->indexToVertexID[data_bar[j]-1] = data_bar[j];
			// //std::cout<<" "<<int(data_bar[j]);
			k++;
		}

		k = 0;
		for (int j=field_offsets[7];j<field_offsets[7]+field_sizes[7];j++){
			b->xVertex[k] = data_bar[field_offsets[8]+k];
			b->yVertex[k] = data_bar[field_offsets[9]+k];
			b->zVertex[k] = data_bar[field_offsets[10]+k];
			k++;
		}

		int l=0;
		for (int j=field_offsets[7];j<field_offsets[7]+field_sizes[7];j++){
			// int m = data_bar[j] - 1;
			for (int k=0; k<3; k++){
				b->cellsOnVertex[l*3+k] = data_bar[field_offsets[12]+3*l+k];
			}
			l++;	
		}


		// MPI_allgather for getting number of cells
		diy::mpi::all_gather(world, nCells, nCell_counts );
		int nCells_global = 0;
		for (auto& n : nCell_counts)
			nCells_global += n;

		diy::mpi::all_gather(world, b->xCell, xCells_global);
		diy::mpi::all_gather(world, b->yCell, yCells_global);
		diy::mpi::all_gather(world, b->zCell, zCells_global);
		mpasobj->xCells.reserve(nCells_global);	
		mpasobj->yCells.reserve(nCells_global);	
		mpasobj->zCells.reserve(nCells_global);	
		for (size_t i=0;i<xCells_global.size();i++){
			mpasobj->xCells.insert(mpasobj->xCells.end(), xCells_global[i].begin(), xCells_global[i].end());

			mpasobj->yCells.insert(mpasobj->yCells.end(), yCells_global[i].begin(), yCells_global[i].end());

			mpasobj->zCells.insert(mpasobj->zCells.end(), zCells_global[i].begin(), zCells_global[i].end());
		}
		mpasobj->nCells = nCells_global;

		fprintf(stderr, "nCell_counts %ld %d %d\n", mpasobj->zCells.size(), nCells_global, nCell_counts[1]);

		// MPI_allgather for getting number of vertices
		diy::mpi::all_gather(world, nVertices, nVertices_counts);
		int nVertices_global = 0;
		for (auto &n : nVertices_counts)
			nVertices_global += n;

		diy::mpi::all_gather(world, b->xVertex, xVertex_global);
		diy::mpi::all_gather(world, b->yVertex, yVertex_global);
		diy::mpi::all_gather(world, b->zVertex, zVertex_global);
		mpasobj->xVertex.reserve(nVertices_global);	
		mpasobj->yVertex.reserve(nVertices_global);	
		mpasobj->zVertex.reserve(nVertices_global);	

		for (size_t i=0;i<xVertex_global.size();i++){
			mpasobj->xVertex.insert(mpasobj->xVertex.end(), xVertex_global[i].begin(), xVertex_global[i].end());

			mpasobj->yVertex.insert(mpasobj->yVertex.end(), yVertex_global[i].begin(), yVertex_global[i].end());

			mpasobj->zVertex.insert(mpasobj->zVertex.end(), zVertex_global[i].begin(), zVertex_global[i].end());
		}
		mpasobj->nVertices = nVertices_global;



		// Construct global kd trees for vertices and cells

		Eigen::MatrixXd C(3, mpasobj->nCells);
		for (int i=0;i<mpasobj->nCells;i++){
			C(0,i) = mpasobj->xCells[i];
			C(1,i) = mpasobj->yCells[i];
			C(2,i) = mpasobj->zCells[i];
		}

		// mpasobj->nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C);

		Eigen::MatrixXd M(3, mpasobj->nVertices);
		for (int i=0;i<mpasobj->nVertices;i++){
			M(0,i) = mpasobj->xVertex[i];
			M(1,i) = mpasobj->yVertex[i];
			M(2,i) = mpasobj->zVertex[i];
		}


		// mpasobj->nns = Nabo::NNSearchD::createKDTreeLinearHeap(M);

		// Construct local kd trees for vertices and cells

		Eigen::MatrixXd C_local(3, nCells);
		for (int i=0;i<nCells;i++){
			C_local(0,i) = b->xCell[i];
			C_local(1,i) = b->yCell[i];
			C_local(2,i) = b->zCell[i];
		}

		mpas_c.nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);

		Eigen::MatrixXd M_local(3, nVertices);
		for (int i=0;i<nVertices;i++){
			M_local(0,i) = b->xVertex[i];
			M_local(1,i) = b->yVertex[i];
			M_local(2,i) = b->zVertex[i];
		}


		mpas_c.nns = Nabo::NNSearchD::createKDTreeLinearHeap(M_local);




	}
	// fprintf(stderr, "indexToCell %f %f\n", data_bar[field_offsets[0]+1], data_bar[field_offsets[0]+0] );

	fprintf(stderr, "xCell %f %f %d\n", b->xCell[0], b->xCell[1], b->indexToCellID[0] );

	fprintf(stderr, "index to cellid size %ld\n", b->indexToCellID.size());
	b->indexToCellID[nCells-1] = 2;

	// populate the recurring members	
	int l=0;
	for (int j=field_offsets[0]; j<field_offsets[0]+field_sizes[0];j++){
		// int m = data_bar[j]-1;
		// fprintf(stderr, " %d", m);
		for (int k=0;k<nVertLevels;k++){
			b->velocityX[l*nVertLevels+k] = data_bar[field_offsets[4]+nVertLevels*l+k];
			b->velocityY[l*nVertLevels+k] = data_bar[field_offsets[5]+nVertLevels*l+k];
			b->velocityZ[l*nVertLevels+k] = data_bar[field_offsets[6]+nVertLevels*l+k];
			b->zTop[l*nVertLevels+k] = data_bar[field_offsets[11]+nVertLevels*l+k];
			b->indexToCellID[l] = data_bar[j];
		}
		l++;
	}


	mpas_c.indexToCellID = b->indexToCellID;
	mpas_c.xCells = b->xCell;
	mpas_c.yCells = b->yCell;
	mpas_c.zCells = b->zCell;
	mpas_c.velocityX = b->velocityX;
	mpas_c.velocityY = b->velocityY;
	mpas_c.velocityZ = b->velocityZ;
	mpas_c.indexToVertexID = b->indexToVertexID;
	mpas_c.xVertex = b->xVertex;
	mpas_c.yVertex = b->yVertex;
	mpas_c.zVertex = b->zVertex;
	mpas_c.zTop = b->zTop;
	mpas_c.cellsOnVertex = b->cellsOnVertex;

	mpas_c.nCells = mpas_c.xCells.size();
	mpas_c.nVertices = mpas_c.xVertex.size();

	mpas_c.nVertLevels = 100; // FIX_ME: Remove the hard coding

	for (int i=0; i<mpas_c.nCells; i++) { // shoudn't i vary from 1 to nCells here? then cellIndex[indexToCellID[i-1]] = i-1
		mpas_c.cellIndex[mpas_c.indexToCellID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
		//	std::cout<<i<<" "<<indexToCellID[i]<<" "<<cellIndex[i]<<"\n";
	}
	for (int i=0; i<mpas_c.nVertices; i++) {
		mpas_c.vertexIndex[mpas_c.indexToVertexID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	}

	// derive velocity on cell verticies (on only one altitude level here)
	mpas_c.velocityXv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.velocityYv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.velocityZv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.zTopVertex.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.zTopVertexNorm.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.xyzCell.resize(mpas_c.nCells*3);
	for (int i=0; i<mpas_c.nCells; i++) mpas_c.xyzCell[i*3] = mpas_c.xCells[i];
	for (int i=0; i<mpas_c.nCells; i++) mpas_c.xyzCell[i*3+1] = mpas_c.yCells[i];
	for (int i=0; i<mpas_c.nCells; i++) mpas_c.xyzCell[i*3+2] = mpas_c.zCells[i];

	for (int i=0; i<mpas_c.nVertices; i++) {
		//    if (cellsOnVertex[i*3] == 0 || cellsOnVertex[i*3+1] == 0 || cellsOnVertex[i*3+2] == 0) continue; // on boundary
		int c0 = mpas_c.cellIndex[mpas_c.cellsOnVertex[i*3]], c1 = mpas_c.cellIndex[mpas_c.cellsOnVertex[i*3+1]], c2 = mpas_c.cellIndex[mpas_c.cellsOnVertex[i*3+2]];

		double X[3][3] = {
			{mpas_c.xyzCell[c0*3], mpas_c.xyzCell[c0*3+1], mpas_c.xyzCell[c0*3+2]},
			{mpas_c.xyzCell[c1*3], mpas_c.xyzCell[c1*3+1], mpas_c.xyzCell[c1*3+2]},
			{mpas_c.xyzCell[c2*3], mpas_c.xyzCell[c2*3+1], mpas_c.xyzCell[c2*3+2]}

		};

		double P[3] = {mpas_c.xVertex[i], mpas_c.yVertex[i], mpas_c.zVertex[i]};

		double lambda[3];
		barycentric_point2triangle(X[0], X[1], X[2], P, lambda);
		//[i * nVertLevels + curVertLevel];
		for (int curVertLevel=0;curVertLevel<mpas_c.nVertLevels; curVertLevel++){

			mpas_c.velocityXv[i*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityX[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityX[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityX[c2*mpas_c.nVertLevels+curVertLevel];
			mpas_c.velocityYv[i*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityY[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityY[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityY[c2*mpas_c.nVertLevels+curVertLevel];
			mpas_c.velocityZv[i*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityZ[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityZ[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityZ[c2*mpas_c.nVertLevels+curVertLevel];

			mpas_c.zTopVertex[mpas_c.nVertLevels*i + curVertLevel] = lambda[0] * mpas_c.zTop[c0*nVertLevels+curVertLevel] + lambda[1] * mpas_c.zTop[c1*nVertLevels+curVertLevel] + lambda[2] * mpas_c.zTop[c2*nVertLevels+curVertLevel];


		}


	}








}

// enqueue remote data
// there is still a link, but you can send to any BlockID = (gid, proc)
void remote_enq(
		PBlock*                              b,
		const diy::Master::ProxyWithLink&   cp,
		const diy::Assigner&                assigner)
{
	// as a test, send my gid to block 1 gids away (in a ring), which is outside the link
	// (the link has only adjacent block gids)
	int my_gid              = cp.gid();
	int dest_gid            = (my_gid + 1) % assigner.nblocks();
	int dest_proc           = assigner.rank(dest_gid);
	diy::BlockID dest_block = {dest_gid, dest_proc};

	// b->values[0] = my_gid;

	// cp.enqueue(dest_block, my_gid);

	std::vector<EndPt> endptvec;
	endptvec.resize(3);

	cp.enqueue(dest_block, endptvec);
}

// dequeue remote data
// there is still a link, but exchange(remote = true) exchanged messages from any block
void remote_deq(PBlock* b, const diy::Master::ProxyWithLink& cp)
{
	std::vector<int> incoming_gids;
	cp.incoming(incoming_gids);
	for (size_t i = 0; i < incoming_gids.size(); i++)
		if (cp.incoming(incoming_gids[i]).size())
		{
			vector<EndPt> recvd_data;
			cp.dequeue(incoming_gids[i], recvd_data);
			fmt::print(stderr, "Remote dequeue: gid {} received value {} from gid {}\n",
					cp.gid(), recvd_data[0].pid, incoming_gids[i]);
			fprintf(stderr, "Received vector: %d\n", 3);
		}
}






// consumer
void con(Decaf* decaf)
{      	
	int max_rounds = 2;
	int max_steps = 1000;

	mpaso mpas_c, mpas_g;
	vector< pConstructData > in_data;


	diy::mpi::communicator world(decaf->con_comm_handle());
	fprintf(stderr, "diy world size: %d", world.size() );

	diy::FileStorage storage("./DIY.XXXXXX");
	int                       nblocks    = world.size();
	int                       threads    = 1;
	int                       mem_blocks = -1;

	diy::Master               master(world,
			threads,
			mem_blocks,
			&create_block,
			&destroy_block,
			&storage,
			&save_block,
			&load_block);


	diy::RoundRobinAssigner assigner(world.size(), nblocks);

	// this example creates a linear chain of blocks with links between adjacent blocks
	std::vector<int> gids;                     // global ids of local blocks
	assigner.local_gids(world.rank(), gids);   // get the gids of local blocks

	PBlock *b;
	int pgid;
	block_io bio;
	for (size_t i = 0; i < gids.size(); ++i)   // for the local blocks in this processor
	{
		int gid = gids[i];

		diy::Link*   link = new diy::Link;   // link is this block's neighborhood
		//        diy::BlockID neighbor;               // one neighbor in the neighborhood
		//        if (gid < nblocks - 1)               // RH neighbors for all but the last block in the global domain
		//        {
		//            neighbor.gid  = gid + 1;                     // gid of the neighbor block
		//            neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
		//            link->add_neighbor(neighbor);                // add the neighbor block to the link
		//        }
		//        if (gid > 0)                         // LH neighbors for all but the first block in the global domain
		//        {
		//            neighbor.gid  = gid - 1;
		//            neighbor.proc = assigner.rank(neighbor.gid);
		//            link->add_neighbor(neighbor);
		//        }

		// keep pointer to current block. use it for populating / updating only current block later
		b = new PBlock;
		pgid = gid;
		global_gid = gid;
		b->gid = gid; // TODO: remove other gid vars e.g. global_gid
		master.add(gid, b, link);    // add the current local block to the master
	}





	std::string ip_file = "graph.topology";
	mpas_g.generate_domain_decomposition_graph(ip_file, world.size());	

	// (remote) exchange some data outside the links
	//	master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
	//		{ 	
	//			remote_enq(b, cp, assigner); }); // seed/deque particles, update vector fields, trace the epoch, and enqueue unfinished particles


	int ctr = 0;
	while (decaf->get(in_data))
	{
		int n_velocityX = 0;
		std::vector<int> indexToCellID;
		std::vector<double> data_bar;
		int fir_cell_idx = 999;
		// get the values and add them
		for (size_t i = 0; i < in_data.size(); i++)
		{

			VectorFliedd d_data_bar = in_data[i]->getFieldData<VectorFliedd>("data_bar");
			if (d_data_bar){
				data_bar = d_data_bar.getVector();
				// printf("got data_bar in con: %f %f %f %f %f %f %f %f %f %f\n", data_bar[13], data_bar[14], data_bar[15], data_bar[16], data_bar[17], data_bar[18], data_bar[1452848+0], data_bar[1452848+4], data_bar[1452848+5], data_bar[1452848+6]);
				// fir_cell_idx = data_bar[0];
				// //printf("fir cell val %f\n", data_bar[int(12+data_bar[4]+data_bar[5]+data_bar[6]+data_bar[7])]);
				// fprintf(stderr, "databar size :%ld\n", data_bar.size());
				// unsigned int init_t=0, fin_t = 40000000, h = 200000;
				// std::vector<double> seeds_xyz(1000*3);
				// std::string op_file = std::to_string(ctr)+".vtp";
				//streamlines slines(mpas1, data_bar, op_file , init_t, fin_t, h, seeds_xyz);

				// use assigner to
				// master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
				//                 {
				//                 populate_block(b, cp, assigner, data_bar ); }); // seed/deque particles, update vector fields, trace the epoch, and enqueue unfinished particles
				int frame_no = (int) data_bar[13];

				if (frame_no==1){

					streamlines slines;
					
					// TODO  premature stopping if work completed before max_rounds 
					populate_block(b, pgid, data_bar, &mpas_g, world, mpas_c);
					int stop = (max_rounds ? max_rounds : 1);
					int incr = (max_rounds ? 1 : 0);
					for (int round = 0; round < stop; round += incr)
					{

						master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
								{

								parallel_streamlines(mpas_g, mpas_c, round, *b, max_steps, cp,  assigner, world);
								});
						bool remote = true;
						master.exchange(remote);	
					}
					bio.write_cell_centers(pgid, world, b->xCell);
					bio.write_particle_traces(pgid, world, *b, max_steps);
				}	

			}else{ 
				fprintf(stderr, "null ptr data_bar\n");
			}


			// foreach load data into block


			// another loop for block synchronous tracing, just keep one traceblock functioni
			// followed by exchange	


			// bool remote = true;
			// master.exchange(remote);
			// master.foreach(&remote_deq);	// particle  tracing here
		}
		// printf("consumer sum = %d\n", sum);
		ctr++;
	}


	// terminate the task (mandatory) by sending a quit message to the rest of the workflow
	fprintf(stderr, "consumer terminating\n");
	decaf->terminate();

}

int main(){


	std::cout<<"starting mpas consumer\n";

	// define the workflow
	Workflow workflow;
	//make_wflow(workflow);
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");

	MPI_Init(NULL, NULL);

	// create decaf
	Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

	// start the task
	con(decaf);

	// cleanup
	delete decaf;
	MPI_Finalize();

	/*
	   std::cout.precision(17);

	   int nseeds = 1;
	//    std::string ip_file = "/Users/mukundraj/Desktop/work/datasets/output.0014-05-01_00.00.00.nc";
	std::string ip_file = "/homes/mraj/work/projects/MPAS-Model-6.0/testdir/MPAS-O_V6.0_QU240/output.nc";
	std::string op_file = "";
	//    unsigned int init_t=0, fin_t = 691199, h = 10000, interval = 172800;
	unsigned int init_t=0, fin_t = 40000000, h = 200000; // iters = fin_t/h
	std::vector<double> seeds_xyz(nseeds*3);

	streamlines slines(ip_file, , init_t, fin_t, h, seeds_xyz);
	//    pathlines plines(ip_file, op_file, init_t, fin_t, interval, h, seeds_xyz);

*/



	std::cout<<"finished\n";

	return 0;

}

bool get_curpos_vel_sl(mpaso &mpas_g, mpaso &mpas_c, double *X, Eigen::Vector3d &c_vel){

	int K = 6;
	double depth = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) - mpas_g.radius;

	// get nearest vertex neighbor ids using mpas_c
	Eigen::VectorXi nearest_idx(K);
	Eigen::VectorXd dists2(K);

	Eigen::VectorXd q(3);
	q<<X[0], X[1], X[2];
	mpas_c.nns->knn(q, nearest_idx, dists2, K);

	// get nearest cell neighbor ids using mpas_c
	Eigen::VectorXi nearest_cell_idx(1);
	Eigen::VectorXd dists2_cell(1);
	mpas_c.nns_cells->knn(q, nearest_cell_idx,dists2_cell, 1);



	// interpolate vertically and get positions and velocity (in values)
	std::vector<double> values;
	values.resize(6*nearest_idx.size());


	double top_boun = mpas_c.zTop[mpas_c.nVertLevels*nearest_cell_idx[0]+0];

	double botom_boun = mpas_c.zTop[mpas_c.nVertLevels*nearest_cell_idx[0]+99];


	if (depth> top_boun || depth< botom_boun){
		// terminate flow line
		return false; // out of global domain
		// fprintf(stderr, "printing values: %f %f %f\n", top_boun, botom_boun, depth);
	}else{ 
		interpolate_vertically(mpas_c.nVertLevels, mpas_c.zTopVertex, nearest_idx, values, depth, mpas_c.xVertex, mpas_c.yVertex, mpas_c.zVertex, mpas_c.velocityXv, mpas_c.velocityYv, mpas_c.velocityZv);

		// interpolate horizontally
		if (!interpolate_horizontally(X[0], X[1], X[2], values, c_vel))
			return false;


	}

	// if (global_gid==1 ){
		// fprintf(stderr,"check vel %f %f %f %f %f %f\n", mpas_c.velocityXv[5539],
		//                                          mpas_c.velocityXv[909],
		//                                           mpas_c.velocityXv[5565],
		//                                            mpas_c.velocityXv[1899],
		//                                             mpas_c.velocityXv[1880],
		//                                              mpas_c.velocityXv[7875]);

	


		//cout<<"inside  get_vel "<<nearest_idx[0]<<" "<<nearest_idx[1]<<" "<<nearest_idx[2]<<" "<<nearest_idx[3]<<" "<<nearest_idx[4]<<" "<<nearest_idx[5]<<" "<<c_vel<<" "<<global_gid<<" "<<mpas_c.velocityXv[5539]<<" "<<mpas_c.velocityXv[909]<<" "<<mpas_c.velocityXv[5565]<<" "<<mpas_c.velocityXv[1899]<<" "<<mpas_c.velocityXv[1880]<<" "<<mpas_c.velocityXv[7875]<<"\n";
	// }
	return true;

}

bool rk1_mpas(mpaso &mpas_g,
		mpaso &mpas_c,
		double *X,
		double h,
		double *Y,
		bool &in_global_domain)
{

	// check if inside global domain
	Eigen::Vector3d c_vel;
	if (get_curpos_vel_sl(mpas_g, mpas_c, X, c_vel)){	

		// advect
		Y[0] = X[0] + h*c_vel[0];	
		Y[1] = X[1] + h*c_vel[1];
		Y[2] = X[2] + h*c_vel[2];
		return true;

	}else{
		in_global_domain = false;
		return false; // needs to be true
	}	
}

void parallel_streamlines(mpaso &mpas_g, mpaso &mpas_c, int round, PBlock &b, int max_steps, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner, diy::mpi::communicator &world){

	// fprintf(stderr, "inside parallel %f\n", mpas_c.xCell_p[0]);
	vector<EndPt> particles;
	map<diy::BlockID, vector<EndPt> > outgoing_endpts;

	double dc_t, dc_c, ratio; //dist from center of earth of top and current layer

	// if first iteration, then seed

	if (round==0){
		int skipval = 5000;
		// std::vector<int> seed_level_ids = {10, 30, 50, 70, 90};
		std::vector<int> seed_level_ids = { 50};
		int n_seed_verts = mpas_c.nVertices/skipval;
		int n_seed_levels = seed_level_ids.size();

		int n_particles = n_seed_verts * n_seed_levels;
		diy::mpi::all_gather(world, n_particles, b.global_nP);
		b.global_start_ids.resize(b.global_nP.size());
		// fprintf(stderr, "global np size %ld \n", b.global_nP.size());
		
		int n_particles_global = b.global_nP[0];
		for (int k = 1; k<b.global_nP.size(); k++){
			b.global_start_ids[k] = b.global_start_ids[k-1]+b.global_nP[k-1];
			// fprintf(stderr, "global_start_pos %d %d\n", k, b.global_start_pos[k]);
			n_particles_global += b.global_nP[k];
		}
		b.global_trace_sizes.resize(n_particles_global);
		fprintf(stderr, "n_seed_verts %d\n", n_seed_verts);
		for (int i=0;i<n_seed_verts; i++){

			for (int j=0;j<n_seed_levels;j++){

				int vert_id = i*skipval;
				if (b.gid==0){
					vert_id = 2230;
				}

				dc_c = mpas_c.radius + mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+seed_level_ids[j]];
				dc_t = mpas_c.radius + mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+0];
				ratio = dc_c/dc_t;
				// fprintf(stderr, "ratio %f\n", ratio);

				EndPt p;
				p.pid = mpas_c.init;
				p.sid = mpas_c.init;
				p.step = 0;
				p.gpid = b.global_start_ids[b.gid] + p.pid;
				p[0] = mpas_c.xVertex[vert_id] * ratio;
				p[1] = mpas_c.yVertex[vert_id] * ratio;
				p[2] = mpas_c.zVertex[vert_id] * ratio;
				if (b.gid==0){
					p[0] = -2089982.455199 * ratio;
					p[1] = 3032836.244409 * ratio;
					p[2] = -5198695.665453 * ratio;

				}


				particles.push_back(p);
				mpas_c.init++;



			}
		}
	
		

	}

	
	// get incoming endpts 	
	vector<int> incoming_gids;
	cp.incoming(incoming_gids);
	for (size_t i = 0; i < incoming_gids.size(); i++)
	{
		if (cp.incoming(incoming_gids[i]).size())
		{
			vector<EndPt> incoming_endpts;
			cp.dequeue(incoming_gids[i], incoming_endpts);
			for (size_t j = 0; j < incoming_endpts.size(); j++) {
				incoming_endpts[j].sid++;
				particles.push_back(incoming_endpts[j]);

			}


		}

	}	


	// trace particles
	int nenq_particles = 0;

	for (int i = 0; i < particles.size(); i++)
	{
		Pt&     cur_p = particles[i].pt; // current end point
		Segment s(particles[i]);         // segment with one point p (initialized with EndPt)
		Pt      next_p;                  // coordinates of next end point
		bool finished = false, in_global_domain = true;
		// int cur_start = s.start_step;
		// int cur_segment_steps = 0;
		// trace particle as far as it will go
		
		while(rk1_mpas(mpas_g, mpas_c, cur_p.coords, 100000, next_p.coords, in_global_domain)){	
			s.pts.push_back(next_p);
			cur_p = next_p;
			if (s.pts.size()%100==0)
				fprintf(stderr, " %d %ld\n", b.gid, s.pts.size()  );
			// cur_segment_steps++;
			// if (s.pts.size() >= max_steps) // CHECK cur_step here instead of s.pts.size()?
			if (s.start_step + s.pts.size() >= max_steps) 
			{
				finished = true;
				break;

			} 


		}
		
		// fprintf(stderr, "pid size %d %ld\n", s.pid, b.trace_sizes.size());
		// b.trace_sizes[s.pid] = s.start_step + s.pts.size();
		b.global_trace_sizes[s.gpid] = s.start_step + s.pts.size();

		// fprintf(stderr, "pts.size %ld\n", s.pts.size());

		if (in_global_domain==false)
			finished = true;


		b.segments.push_back(s);	

		if (finished){
			mpas_c.done++;

			fprintf(stderr,"finished particle\n");
		}
		else{ 			// package unfinished endpoint for sending 
			fprintf(stderr,"enqueuing particle\n");
			EndPt out_pt(s); // TO_DO: Also update the start_time here
			int dest_gid = mpas_c.get_bid_for_pt(next_p.coords);		
			int dest_proc = assigner.rank(dest_gid);

			diy::BlockID bid = {dest_gid, dest_proc};

			outgoing_endpts[bid].push_back(out_pt);


			nenq_particles++;

		}	

	}	
	

	// send off unfinished segments using rexchange

	for (map<diy::BlockID, vector<EndPt> >::const_iterator it =
			outgoing_endpts.begin(); it != outgoing_endpts.end(); it++)
		cp.enqueue(it->first, it->second);

}


void block_io::write_particle_traces(int gid, const diy::mpi::communicator &world, PBlock &b, int max_steps){

	int ret, ncfile, nprocs, rank;
	char filename[256], buf[13] = "Hello World\n";
	int data[2];


	// int nP = b.trace_sizes.size(); // number of particles in current block
	// std::vector<std::vector<int>> trace_sizes_global;
	// concatenate trace_from all blocs in proc that processes block 0

	// diy::mpi::gather(world, b.trace_sizes, trace_sizes_global, 0);

	// if (rank == 0){
	//         diy::mpi::reduce(world, b.global_trace_sizes, b.global_trace_sizes,0, MPI_MAX);
	// }
	// else{
	//         diy::mpi::reduce(world, b.global_trace_sizes);
	// }
	//
	// // fprintf(stderr, "np %d\n", nP);
	//
	//

	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &nprocs);
	std::vector<int> trace_sizes_recv;
	trace_sizes_recv.resize(b.global_trace_sizes.size());
	// fprintf(stderr, "trace_sizes_global size %ld %d\n", trace_sizes_global.size(), rank);

	int trace_sizes_global[b.global_trace_sizes.size()];
	for (int i=0; i<b.global_trace_sizes.size();i++){
		trace_sizes_global[i] = b.global_trace_sizes[i];
	}	


	if (rank==0){	
		MPI_Reduce(MPI_IN_PLACE, &trace_sizes_global, b.global_trace_sizes.size(), MPI_INT, MPI_MAX, 0, world);
		// for (int i=0; i<b.global_trace_sizes.size();i++){
		//         fprintf(stderr, " %d",trace_sizes_global[i]);
		// }
		// fprintf(stderr, " \n");
	}else{

		MPI_Reduce(&trace_sizes_global, &trace_sizes_global, b.global_trace_sizes.size(), MPI_INT, MPI_MAX, 0, world);
	}


	strcpy(filename, "particle_traces.nc");


	int dimid_s, dimid_p; // varid for steps and particles
	int ndims_pos = 2;
	int varid_xPos, varid_yPos, varid_zPos, varid_trace_sizes_global;
	int pos_dims[ndims_pos];

	ret = ncmpi_create(world, filename,
			NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	// define dimensions
	ret = ncmpi_def_dim(ncfile, "nStep", NC_UNLIMITED, &dimid_s);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_dim(ncfile, "nParticle", b.global_trace_sizes.size(), &dimid_p);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	pos_dims[0] = dimid_s;
	pos_dims[1] = dimid_p;
	// define variables
	ret = ncmpi_def_var(ncfile, "xPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_xPos);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, "yPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_yPos);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, "zPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_zPos);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// MPI_Offset start[ndims_pos], count[ndims_pos];

	// end define mode
	ret = ncmpi_enddef(ncfile);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// write data to file each segment as a separate call

	// TODO create and allocate buffer 
	// http://cucis.ece.northwestern.edu/projects/PnetCDF/doc/pnetcdf-c/ncmpi_005fput_005fvarn_005f_003ctype_003e.html
	MPI_Offset **starts, **counts;

	int num_reqs = b.segments.size();
	/* allocate starts and counts */
	starts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
	starts[0] = (MPI_Offset*)  calloc(num_reqs*ndims_pos, sizeof(MPI_Offset));
	for (int i=1; i<num_reqs; i++)
		starts[i] = starts[i-1] + ndims_pos;

	counts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
	counts[0] = (MPI_Offset*)  calloc(num_reqs*ndims_pos, sizeof(MPI_Offset));
	for (int i=1; i<num_reqs; i++)
		counts[i] = counts[i-1] + ndims_pos;
	
	for (int i=0; i<b.segments.size(); i++){
		starts[0][0] = 0; starts[0][1] = b.segments[i].gpid;
		counts[0][0] = b.segments[i].pts.size(); counts[0][1] = 1;
	}

       	/* allocate write buffer */
	int buf_len = 0;
	for (int i=0; i<num_reqs; i++) {
		MPI_Offset w_req_len=1;
		for (int j=0; j<ndims_pos; j++)
			w_req_len *= counts[i][j];
		buf_len += w_req_len;
	}
	fprintf(stderr, "buflen num_reqs, ndims_pos %d %d %d\n", buf_len, num_reqs, ndims_pos);
	double *buffer = (double*) malloc(buf_len * sizeof(double));
	double *buffer_y = (double*) malloc(buf_len * sizeof(double));
	double *buffer_z = (double*) malloc(buf_len * sizeof(double));
 
	// for (int i=0; i<buf_len; i++) {
	//         buffer[i] = (double)rank;
	//         fprintf(stderr, "bufval %f\n", buffer[i]);
	// }
	
	int bidx = 0;
	for (int i=0; i<b.segments.size(); i++) {
		for (int j=0; j<b.segments[i].pts.size();j++){
			buffer[bidx] = b.segments[i].pts[j].coords[0];
			buffer_y[bidx] = b.segments[i].pts[j].coords[1];
			buffer_z[bidx] = b.segments[i].pts[j].coords[2];
			bidx++;
		}
	}        

	fprintf(stderr, "after allocation\n"); 		

	ret = ncmpi_put_varn_double_all(ncfile, varid_xPos, num_reqs ,starts, counts, buffer);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_put_varn_double_all(ncfile, varid_yPos, num_reqs ,starts, counts, buffer_y);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_put_varn_double_all(ncfile, varid_zPos, num_reqs ,starts, counts, buffer_z);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	for (int i=0;i<b.global_trace_sizes.size();i++){
		fprintf(stderr, "global %d %d\n", trace_sizes_global[i], b.gid);	
	}
	
	// close file
	ret = ncmpi_close(ncfile);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	free(buffer);
	free(buffer_y);
	free(buffer_z);
	// TODO free starts and counts

	
	
	if (rank == 0){

		ret = ncmpi_open(MPI_COMM_SELF, filename,
				NC_WRITE|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		ret = ncmpi_redef(ncfile); /* enter define mode */
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		//	Next define trace_sizes global
		ret = ncmpi_def_var(ncfile, "trace_sizes_global", NC_INT, 1, &dimid_p, &varid_trace_sizes_global);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

	
		// end define mode
		ret = ncmpi_enddef(ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		fprintf(stderr, "here \n");

		MPI_Offset start[1], count[1];
		start[0]=0, count[0]=b.global_trace_sizes.size();
		ret = ncmpi_put_vara_int_all(ncfile, varid_trace_sizes_global, start, count, trace_sizes_global);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);	
		fprintf(stderr, "crossed \n");

		// close file
		ret = ncmpi_close(ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);


	}
}
