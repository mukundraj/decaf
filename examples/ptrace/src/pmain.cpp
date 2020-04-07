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

#include "utils/opts.h"

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
			// dprint("incoming in %d", cp.gid());
        }
       
    }
   
}


bool trace_particles(block *b,
                     const diy::Master::ProxyWithLink &cp,
                     const diy::Assigner &assigner, 
					 const int max_steps, 
					 pathline &pl, 
					 int prediction, 
					 size_t &nsteps, 
					 std::vector<EndPt>& particles_hold){
				bool val = true;

				// b->particles_store.clear();
				deq_incoming_iexchange(b, cp);
				// dprint("calling trace_particles");
				val = pl.compute_streamlines(b, cp, assigner, prediction, nsteps, particles_hold);

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
	int check = 0;                  // write out traces
	int gid;

	int ndims = 3;           // domain dimensions
	string particle_file;    // input file name
	string gp_file;    		// input file name
	double dtSim = 5000*7200*10;
	double dtParticle = 300000;
	size_t nCells;

	bool prediction = true;
	double time_prep=0, time_predrun=0, time_kdtree=0, time_readdata=0, time_filter=0, time_final=0, time_predrun_loc=0, time_final_loc, time_trace = 0;


	using namespace opts;

	 // command-line ags
    Options ops(argc, argv);

	ops >> Option('b', "blocks", nblocks, "Total number of blocks to use")
		>> Option('c', "check", check, "Write out traces for checking");

	 if (ops >> Present('h', "help", "show help") ||
        !(ops >> PosOption(particle_file) >> PosOption(gp_file) >> PosOption(prediction) >> PosOption(dtSim) >> PosOption(dtParticle)))
    {
        if (world.rank() == 0)
        {
            fprintf(stderr, "Check ops usage \n", argv[0]);
            cout << ops;
        }
        return 1;
    }

	size_t nsteps = 0, ntransfers = 0, nsteps_lagged = 0, nsteps_pred=0;

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
		std::string fname_graph = "graph.info", fname_graphpart = gp_file;
		//"graph.info.part." + std::to_string(world.size());
		dprint("fname_graphpart %s", fname_graphpart.c_str());

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
		std::string fname_particles = particle_file; //"particles.nc";
		b->init_seeds_particles(world, fname_particles, 0);

		b->init_partitions(); // must be called after loading mpas data, after setting b->gid, and after init_seeds_partitions
		nCells = b->gcIdxToGid.size();
		

		master.add(gid, b, link);
	}
	// prediction advection using iexchange
	// double dtSim = 7200, dtParticle = 300;
	// double dtSim = 5000*7200*10, dtParticle = 300000;

	
	std::vector<EndPt> particles_hold; // staging area for holding to-be-predicted particles

	if (prediction)
	{	
		Bounds domain{DIM};
		for (unsigned i = 0; i < DIM; ++i)
		{
			domain.min[i] = -6371329. - 100; // 6371229.
			domain.max[i] = 6371329. +100;
		}

		
		diy::Master master_kdt(world,
					   threads,
					   mem_blocks,
					   NULL,
					   NULL,
					   &storage,
					   NULL,
					   NULL);
	// 		std::vector<int> gids;                     // global ids of local blocks
    // assigner.local_gids(world.rank(), gids);   // get the gids of local blocks

		for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
		{	block* b = new block();
			int gid = gids[i];
			b->gid = gid;

			RCLink*         l   = new RCLink(DIM, domain, domain);
			master_kdt.add(gid, b, l);

		}

		

		world.barrier();
        double time0 = MPI_Wtime(); 
		
		// select prediction particles; store tag along particles in b->particles_hold
		master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			std::random_device rd;
			std::mt19937 g(4);
			std::shuffle(b->particles.begin(), b->particles.end(), g);

			size_t pred_size = b->particles.size()/1;

			// b->particles_hold.insert(std::end(b->particles_hold), std::begin(b->particles) + pred_size, std::end(b->particles));
			particles_hold.insert(std::end(particles_hold), std::begin(b->particles) + 0, std::end(b->particles));
			b->particles.resize(pred_size);

			for (EndPt& ep: particles_hold){ // updating p{xyz} to cell{xyz} for load based sorting
				ep.pt_hold = ep.pt;
				ep[0] = b->xCell[ep.glCellIdx];
				ep[1] = b->yCell[ep.glCellIdx];
				ep[2] = b->zCell[ep.glCellIdx];
			}

		});

		world.barrier();
		double time1 = MPI_Wtime();
		time_prep = time1 - time0;

		// prediction advection that appends work locations to particles_hold

		// master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {

		


		// });

		// prediction advection using iexchange
		master.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {

			pathline pl(*b, dtSim, dtParticle);

			// update_velocity_vectors: both timesteps point to same vectors
			pl.set_velocity_vectors(*b);


			bool val = trace_particles(b,
									   icp,
									   assigner,
									   max_steps,
									   pl,
									   true,
									   nsteps,
									   particles_hold);

			return val;
		});

		block* block_ptr;
		// replace particle p{xyz} with cell{xyz}, move b->particles_hold to b->particles
		master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {

			for (auto i: b->in_partition_set){
				EndPt cell_cen;
				cell_cen[0] = b->xCell[i];
				cell_cen[1] = b->yCell[i];
				cell_cen[2] = b->zCell[i];
				cell_cen.predonly = 2; // indicates this is the cell center point
				cell_cen.glCellIdx = i;
				cell_cen.source_gid = b->gid;
				particles_hold.push_back(cell_cen);
				if (cell_cen.glCellIdx==196056 || cell_cen.glCellIdx==48616){
					dprint("cell_cen.glCellIdx %d starts from gid %d", cell_cen.glCellIdx , b->gid);
				}
			}

			// b->particles = std::move(b->particles_hold);
			// dprint("particles size %ld", b->particles.size());
			block_ptr = std::move(b);
		});

		master_kdt.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			b->particles = std::move(particles_hold);
		});
		
		// load based repartition
		bool wrap = false;
		size_t samples = 512;
		diy::kdtree(master_kdt, assigner, DIM, domain, &block::particles, samples, wrap);

		master_kdt.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			particles_hold = std::move(b->particles);
			gid = b->gid;
		});

		dprint("Populating gcIdxToGid_global. nCells %d", nCells);
		// preparing to compute links by computing gcIdxToGid_global	
		std::vector<int> gcIdxToGid_local(nCells);
		std::vector<int> gcIdxToGid_global(nCells);

		for (size_t i=0; i<particles_hold.size(); i++){
				if (particles_hold[i].predonly==2){
					gcIdxToGid_local[particles_hold[i].glCellIdx] = gid;
				}
		}
		diy::mpi::all_reduce (world, gcIdxToGid_local, gcIdxToGid_global,  std::plus<int>());


		dprint("Finished prediction .............. ");
		


		// post prediction run with new master with reconstructed links
		nsteps = 0;

		diy::Master master_final(world,
							   threads,
							   mem_blocks,
							   NULL,
							   NULL,
							   &storage,
							   NULL,
							   NULL);

		for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
		{
			block *b = std::move(block_ptr);
			int gid = gids[i];
			b->gid = gid;

			b->gcIdxToGid = std::move(gcIdxToGid_global);
			std::string fname_graph = "graph.info";
			std::set<int> new_neighbors;
			b->create_links_from_gcIdxToGid(fname_graph, new_neighbors);

			diy::Link *link = new diy::Link; // link is this block's neighborhood

			// dprint("new_neighbors %ld", new_neighbors.size());
			for (int lgid : new_neighbors)
			{
				// dprint("rank %d links to %ld", world.rank(), i);
				diy::BlockID neighbor;
				neighbor.gid = lgid;						 // gid of the neighbor block
				neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
				link->add_neighbor(neighbor);
				// dprint("gid %d nbr %d", gid, lgid);
			}
			master_final.add(gid, b, link);
		} 



		master_final.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			 // read mpas data
			std::string fname_data = "output.nc";
			// b->loadMeshFromNetCDF_CANGA(world, fname_data, 0);

			// // repeating to populate gcIdxToGid, don't use the links identified here
			// std::string fname_graph = "graph.info", fname_graphpart = "graph.info.part." + std::to_string(world.size());
			// set<int> links;
			// b->create_links(fname_graph, fname_graphpart, links);

		});

		dprint("Read data again done ...... " );

		// filter out the cell center particles and workload particles, update b->in_partition


		
		master_final.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {

			b->particles = std::move(particles_hold);
			b->in_partition.resize(b->gcIdxToGid.size());
			
			// std::fill(b->in_partition.begin(), b->in_partition.end(), 0);
			// std::fill(b->gcIdxToGid.begin(), b->gcIdxToGid.end(), 0);
			
			

			for (size_t i=0; i<b->particles.size(); i++){
				if (b->particles[i].predonly==2){
					b->in_partition[b->particles[i].glCellIdx] = 1;
					if (b->particles[i].glCellIdx==196056 || b->particles[i].glCellIdx==48616 ){
						dprint("cell_cen.glCellIdx %d goes to gid %d", b->particles[i].glCellIdx , b->gid);
					}
				}else if (b->particles[i].predonly==0){
					b->particles[i].pt = b->particles[i].pt_hold;
					particles_hold.push_back(b->particles[i]);
				}
			}

			

			// dprint("post partition %d, particles %ld", postsum, b->particles.size());
			b->particles = std::move(particles_hold);
		});

		dprint("Starting second advection .............. ");


			

			// final advection using iexchange
			master_final.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {
				
				pathline pl(*b, dtSim, dtParticle);

				// update_velocity_vectors: both timesteps point to same vectors
				pl.set_velocity_vectors(*b);

				bool val = trace_particles(b,
											icp,
											assigner,
											max_steps,
											pl, 
											false,
											nsteps,
											particles_hold);

				

				
				return val;
			});

			dprint("done final advection");
			


	}else{

		// prediction advection using iexchange
		master.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {
			pathline pl(*b, dtSim, dtParticle);

			// update_velocity_vectors: both timesteps point to same vectors
			pl.set_velocity_vectors(*b);

			bool val = trace_particles(b,
									   icp,
									   assigner,
									   max_steps,
									   pl,
									   false,
									   nsteps,
									   particles_hold);

			return val;
		});

		// dprint("rank %d, segs %ld", world.rank(), b->segments.size());

		master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			if (b->segments.size() == 0)
			{
				Pt p;
				p.coords[0] = 0;
				p.coords[1] = 0;
				p.coords[2] = 0;

				Segment seg;
				seg.pid = -1;
				seg.pts.push_back(p);
				b->segments.push_back(seg); // just dealing with the quirk in pnetcdf writer if empty block happens to be the last block
			}
		});
	}
    


	if (check==1 && prediction == 0)
    {
		// write out segments
		master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			b->parallel_write_segments(world, 0);
		});
	}


	// write out stats
	size_t nsteps_global=0;
	diy::mpi::reduce(world, nsteps, nsteps_global, 0, std::plus<size_t>());
	size_t maxsteps_global=0;
    diy::mpi::reduce(world, nsteps, maxsteps_global, 0, diy::mpi::maximum<size_t>());

	float avg = float(nsteps_global) / world.size();
    float balance = (float(maxsteps_global))/ float(avg);

	double  time_predrun_loc_max, time_predrun_loc_avg, time_fin_loc_max, time_fin_loc_avg, time_trace_max, time_trace_avg;
	size_t minsteps_global, ntransfers_global;
	if (world.rank() == 0)
	{
		fprintf(stderr, "predd , %d, nsteps_global , %ld, maxsteps_global , %ld, bal , %f, time_tot , %f, time_overhead, %f, worldsize, %d, minsteps, %ld,\n", prediction, nsteps_global, maxsteps_global, balance, double(0.0), double(0.0), world.size(), double(0.0));
		dprint("times: predrun, %f, kdtree , %f, readdata, %f, filter ,%f, final , %f, prediction, %d, max, %ld, min, %ld, nsteps, %ld, wsize, %d, time_pred, %f, tot_transfers, %ld, prdrun_local(max avg), %f, %f, fin_local (max avg), %f, %f, max_steps, %d, time_trace_max, %f, time_trace_avg, %f, ", time_predrun, time_kdtree, time_readdata, time_filter, time_final, prediction, maxsteps_global, minsteps_global, nsteps_global, world.size(), time_prep, ntransfers_global, time_predrun_loc_max, time_predrun_loc_avg, time_fin_loc_max, time_fin_loc_avg, max_steps, time_trace_max, time_trace_avg);
		dprint("done");
	}


   
}
