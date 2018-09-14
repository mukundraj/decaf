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

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/serialization.hpp>


void populate_block( PBlock* b,
		int gid,
		std::vector<double> &data_bar ){


	int frame_no = (int)data_bar[13];
	int nVertices = (int) data_bar[14];
	int nCells = (int) data_bar[15];
	int nVertLevels = (int) data_bar[16];

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
		b->indexToCellID.resize(nCells); //0
		b->xCell.resize(nCells); //1
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

		b->xCell[0] = 2.9;
		b->xCell[nCells-1] = 9.2;
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


int frame_number = 0;




// consumer
void con(Decaf* decaf)
{       
	mpaso mpas1;
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
		master.add(gid, b, link);    // add the current local block to the master
	}

	int num_iters = 2;
	for (int i = 0; i < num_iters; i++)
	{


	}


	std::string ip_file = "graph.topology";
	mpas1.generate_domain_decomposition_graph(ip_file, world.size());	

	// (remote) exchange some data outside the links
	master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
			{ 	
			remote_enq(b, cp, assigner); }); // seed/deque particles, update vector fields, trace the epoch, and enqueue unfinished particles


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
				populate_block(b, pgid, data_bar);


			}else{ 
				fprintf(stderr, "null ptr data_bar\n");
			}


			// foreach load data into block


			// another loop for block synchronous tracing, just keep one traceblock functioni
			// followed by exchange	


			bool remote = true;
			master.exchange(remote);
			master.foreach(&remote_deq);	// particle  tracing here
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
