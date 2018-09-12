
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
#include "ptrace.hpp"

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/serialization.hpp>

struct TBlock
{
    TBlock(): count(0)                   {}

    int   count;
    vector<int> values{9,9,9};

};

namespace diy
{
    template<>
        struct Serialization<TBlock>
	{
        static void save(BinaryBuffer& bb, const TBlock& b)
	{
            diy::save(bb, b.count);
            diy::save(bb, b.values);
        
	}

        static void load(BinaryBuffer& bb, TBlock& b)
	{
            diy::load(bb, b.count);
            diy::load(bb, b.values);
        
	}
    
	};

}

void* create_block()                      { return new TBlock;  }
void  destroy_block(void* b)              { delete static_cast<TBlock*>(b);  }
void  save_block(const void* b,
                 diy::BinaryBuffer& bb)   { 
						diy::save(bb, *static_cast<const TBlock*>(b));
											
						  
					}
void  load_block(void* b,
                 diy::BinaryBuffer& bb)   { diy::load(bb, *static_cast<TBlock*>(b));  }


// enqueue remote data
// there is still a link, but you can send to any BlockID = (gid, proc)
void remote_enq(
        TBlock*                              b,
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
void remote_deq(TBlock* b, const diy::Master::ProxyWithLink& cp)
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

        master.add(gid, new TBlock, link);    // add the current local block to the master
    }

       int num_iters = 2;
    for (int i = 0; i < num_iters; i++)
    {

        
    }

	
	std::string ip_file = "graph.topology";
	mpas1.generate_domain_decomposition_graph(ip_file, world.size());	

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
				// data_bar = d_data_bar.getVector();
				// printf("got data_bar in con: %f %f %f %f %f %f %f %f %f %f\n", data_bar[13], data_bar[14], data_bar[15], data_bar[16], data_bar[17], data_bar[18], data_bar[1452848+0], data_bar[1452848+4], data_bar[1452848+5], data_bar[1452848+6]);
				// fir_cell_idx = data_bar[0];
				// //printf("fir cell val %f\n", data_bar[int(12+data_bar[4]+data_bar[5]+data_bar[6]+data_bar[7])]);
				// fprintf(stderr, "databar size :%ld\n", data_bar.size());
				// unsigned int init_t=0, fin_t = 40000000, h = 200000;
				// std::vector<double> seeds_xyz(1000*3);
				// std::string op_file = std::to_string(ctr)+".vtp";
				//streamlines slines(mpas1, data_bar, op_file , init_t, fin_t, h, seeds_xyz);
				
				// use assigner to 

			}else{ 
				fprintf(stderr, "null ptr data_bar\n");
			}


			

		// (remote) exchange some data outside the links
		master.foreach([&](TBlock* b, const diy::Master::ProxyWithLink& cp)
			{ 	
				remote_enq(b, cp, assigner); }); // seed/deque particles, update vector fields, trace the epoch, and enqueue unfinished particles
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
