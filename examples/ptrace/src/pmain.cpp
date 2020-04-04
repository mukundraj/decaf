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

// #include "misc.h"

#include "mpas_io.h"
#include "pathline.h"
#include "block.h"
#include "misc.h"
#include <diy/algorithms.hpp>
#include <cmath>

typedef diy::RegularContinuousLink RCLink;
typedef diy::ContinuousBounds Bounds;
static const unsigned DIM = 3;



using namespace std;


void deq_incoming_iexchange(block *b,
                            const diy::Master::ProxyWithLink &cp)
{
    diy::Link *l = static_cast<diy::Link *>(cp.link());
    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        while (cp.incoming(nbr_gid)){
            EndPt incoming_endpt;
            cp.dequeue(nbr_gid, incoming_endpt);
            b->particles_store.push_back(incoming_endpt);
			dprint("incoming in %d", cp.gid());
        }
       
    }
   
}


bool trace_particles(block *b,
                     const diy::Master::ProxyWithLink &cp,
                     const diy::Assigner &assigner, 
					 const int max_steps, 
					 pathline &pl){
				bool val = true;

				// b->particles_store.clear();
				deq_incoming_iexchange(b, cp);
				// dprint("calling trace_particles");
				val = pl.compute_streamlines(b, cp, assigner);

				return val;

}

int main(int argc, char* argv[])
{
     diy::mpi::environment     env(argc, argv); // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;

    diy::FileStorage storage("./DIY.XXXXXX");
	int nblocks = world.size();
	int threads = 1;
	int mem_blocks = -1;
	int max_steps = 5;

	diy::Master master(world,
					   threads,
					   mem_blocks,
					   NULL,
					   NULL,
					   &storage,
					   NULL,
					   NULL);

	diy::RoundRobinAssigner assigner(world.size(), nblocks);


	// add master to block 
	std::vector<int> gids;                     // global ids of local blocks
    assigner.local_gids(world.rank(), gids);   // get the gids of local blocks

	for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
    {	block* b = new block();
		int gid = gids[i];
		b->gid = gid;

		diy::Link*    link = new diy::Link;  // link is this block's neighborhood


		// read mpas data
		std::string fname_data = "output.nc";
		b->loadMeshFromNetCDF_CANGA(world, fname_data, 0);
		

		
		// read and add link
		std::string fname_graph = "graph.info", fname_graphpart = "graph.info.part." + std::to_string(world.size());
        set<int> links;
        b->create_links(fname_graph, fname_graphpart, links);
		for (int lgid: links){
			// dprint("rank %d links to %ld", world.rank(), i);
			diy::BlockID  neighbor;
			neighbor.gid  = lgid;                     // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			link->add_neighbor(neighbor);
		}
		// init particles
		std::string fname_particles = "particles.nc";
		b->init_seeds_particles(world, fname_particles, 0);

		b->init_partitions(); // must be called after loading mpas data, after setting b->gid, and after init_seeds_partitions
		

		master.add(gid, b, link);
	}
	// prediction advection using iexchange
	// double dtSim = 7200, dtParticle = 300;
	double dtSim = 5000*7200*10, dtParticle = 300000;
	master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {

		pathline pl(*b, dtSim, dtParticle);

		// update_velocity_vectors: both timesteps point to same vectors
		pl.set_velocity_vectors(*b);


		// prediction advection using iexchange
		master.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {
			
			bool val = trace_particles(b,
										icp,
										assigner,
										max_steps,
										pl);

			
			return val;
		});

		if (b->segments.size()==0){
			Pt p;
			p.coords[0] = 0; 
			p.coords[1] = 0; 
			p.coords[2] = 0; 

			Segment seg;
			seg.pid = -1;
			seg.pts.push_back(p);
			b->segments.push_back(seg); // just dealing with the quirk in pnetcdf writer if empty block happens to be the last block
		}

		dprint("rank %d, segs %ld", world.rank(), b->segments.size());
			
		

	});

	// rebalance
	
    
    // final advect using iexchange


    // write out segments
	master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
		b->parallel_write_segments(world, 0);
	});

    dprint("done");
}
