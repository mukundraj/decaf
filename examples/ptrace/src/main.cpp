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
#include "core/pathline.h"
#include "core/streamline.h"
#include "utils/process.h"

typedef diy::RegularContinuousLink RCLink;
typedef diy::ContinuousBounds Bounds;
static const unsigned DIM = 3;








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
					 flowline *fl, 
					 int prediction, 
					 size_t &nsteps, 
					 std::vector<EndPt>& particles_hold, 
					 int skip_rate){
				bool val = true;

				// b->particles_store.clear();
				deq_incoming_iexchange(b, cp);
				// dprint("calling trace_particles");
				val = fl->compute_flow(b, cp, assigner, prediction, nsteps, particles_hold, skip_rate);

				return val;

}



bool first_done = false;

size_t nsteps = 0;
double time_prep=0, time_predrun=0, time_kdtree=0, time_readdata=0, time_filter=0, time_final=0;

// consumer
void con(Decaf *decaf, block &mpas1, pathline &pl, string &gp_file, double dtSim, double dtParticle, std::string& particle_file, int seed_rate, int pred_percent, const int max_steps, int skip_rate, int prediction, int nblocks, int rebal_interval)
{
	std::string dbgmsg;

	// mpaso mpas_c;
	vector<pConstructData> in_data;
	diy::mpi::communicator world(decaf->con_comm_handle());

	// recomputed particle file ("particles2.nc") includes global cell indices for seeds
	std::string fname_particles = "particles.nc", fname_graphinfo = "graph.info.part." + std::to_string(world.size());


	diy::FileStorage storage("./DIY.XXXXXX");
	int threads = 2;
	int mem_blocks = -1;

	Bounds domain{DIM};
	for (unsigned i = 0; i < DIM; ++i)
	{
		domain.min[i] = -6371329.; // 6371229.
		domain.max[i] = 6371329.;
	}

	block * block_ptr = NULL;
	set<int> links;
	std::vector<EndPt> particles_hold;
	diy::Link*    link_pred = new diy::Link;  // link is this block's neighborhood

	// diy::Master master_pred(world,
	// 				   threads,
	// 				   mem_blocks,
	// 				   NULL,
	// 				   NULL,
	// 				   &storage,
	// 				   NULL,
	// 				   NULL);

	diy::RoundRobinAssigner assigner(world.size(), nblocks);

	std::vector<int> gids;					 // global ids of local blocks
	assigner.local_gids(world.rank(), gids); // get the gids of local blocks

	std::string fname_graph = "graph.info", fname_graphpart = gp_file;
	std::string fname_data = "output2.nc";

	for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
    {	
		block* b = new block();
		int gid = gids[i];
		b->gid = gid;
		b->loadMeshFromNetCDF_CANGA(world, fname_data, 0, fname_graph);

	
		b->init_seeds_particles(world, fname_particles, 0, seed_rate);	
		dprint("particles %ld, rank %d", b->particles.size(), gid);

		// b->init_partitions(); // must be called after loading mpas data, after setting b->gid, and after init_seeds_partitions; populates b->in_partition with bools and in_partition_set

		// read and add link
        set<int> links;
        b->create_links(fname_graph, fname_graphpart, links);

		b->gcIdToGid_local.resize(b->nCells + 1);
		b->gcIdToGid.resize(b->nCells + 1);
		b->gcIdToGid_init.resize(b->nCells + 1);

		
		

		for (int lgid: links){
			// dprint("rank %d links to %ld", world.rank(), i);
			diy::BlockID  neighbor;
			neighbor.gid  = lgid;                     // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			link_pred->add_neighbor(neighbor);
		}

		// master_pred.add(block_ptr->gid, block_ptr, link_pred);

		block_ptr = std::move(b);

	}

	// creating the master rebalance and advection
	diy::Master master_baladv(world,
					   threads,
					   mem_blocks,
					   NULL,
					   NULL,
					   &storage,
					   NULL,
					   NULL);

	RCLink *l = new RCLink(DIM, domain, domain);
	master_baladv.add(block_ptr->gid, block_ptr, l);

	dbgmsg = "###particles outside loop";
	print_global_num_particles(dbgmsg, world, master_baladv, 0);



	while (decaf->get(in_data))
	{	
		std::vector<double> data_bar_dbl;
		std::vector<int> data_bar_int;
		int data_id, framenum;

		// get the values and add them
		for (size_t i = 0; i < in_data.size(); i++)
		{

			SimpleFieldi d_metadata = in_data[i]->getFieldData<SimpleFieldi>("field_id");
			// SimpleFieldi d_metadata = in_data.at("in")->getFieldData<SimpleFieldi>("field_id");
			SimpleFieldi frame_num = in_data[i]->getFieldData<SimpleFieldi>("frame_num");
			data_bar_int = in_data[i]->getFieldData<VectorFieldi>("data_i").getVector();
			data_bar_dbl = in_data[i]->getFieldData<VectorFliedd>("data_d").getVector();

			
			data_id = d_metadata.getData();
			framenum = frame_num.getData();
			if (world.rank()==0){
				dprint("data_id %d, framenum %ld, %ld %ld, %d %f", data_id, framenum, data_bar_int.size(), data_bar_dbl.size(), data_bar_int[0], data_bar_dbl[0]);
			}

			block_ptr->update_data(data_id, framenum, data_bar_int, data_bar_dbl);
			

			
		}

		if ((data_id == 17 && framenum > 1) || (data_id == 15 && framenum == 1))
		{
			if (world.rank()==0)
				dprint("starting new frame , %d", framenum);

			if ((framenum - 1) % rebal_interval == 0){

				if (world.rank()==0)
					dprint("starting prediction in frame , %d", framenum);

				particles_hold.clear();

				block_ptr->position_seeds_before_prediction(pred_percent, particles_hold);
				// dprint("block_ptr->particles %d", block_ptr->particles.size());
				block_ptr->position_data_before_prediction(world, block_ptr->gid, framenum, particles_hold);

				// dprint("particles_hold %d", particles_hold.size());
				// for(size_t i=0; i<particles_hold.size(); i++){
				// 	if (particles_hold[i].predonly==2)
				// 		dprint("found 2");
				// 	// else dprint("not 2 but %d", particles_hold[i].predonly);

				// }

				// if (world.rank()==0){
				// 	dprint("gcIdToGid_init %ld", block_ptr->gcIdToGid_init.size());
				// 	pvi(block_ptr->gcIdToGid_init);
				// }

				
				
				
				// ghost exchange
				diy::Master master_pred(world,
					   threads,
					   mem_blocks,
					   NULL,
					   NULL,
					   &storage,
					   NULL,
					   NULL);
					   

				std::vector<int> gids;					 // global ids of local blocks
				assigner.local_gids(world.rank(), gids); // get the gids of local blocks
				diy::Link*    link = new diy::Link;  // link is this block's neighborhood
				block_ptr->update_link(links, assigner, link);				
				master_pred.add(gids[0], block_ptr, link);

				
				// if (world.rank()==0)

				// ghost exchange
				ghost_exchange(master_pred, assigner, true);

				int ctr = 2630-1;
				if (block_ptr->gid==0)
					for (size_t i=ctr*block_ptr->nVertLevels; i<ctr*block_ptr->nVertLevels +block_ptr->nVertLevels; i++)
						fprintf(stderr, "%f ", block_ptr->velocitiesX[1][i]);


				// if (world.rank()==0)
				// dprint("After ghost exchange, local_gcIds %ld, framenum %d", block_ptr->local_gcIds.size(), framenum);		

				world.barrier();

				// interpolate to vertex
				interpolate_cell_to_vertex(master_pred, true);
				world.barrier();

				ctr = 869;
				if (block_ptr->gid==0){
					dprint(" ");
					for (size_t i=ctr*block_ptr->nVertLevels; i<ctr*block_ptr->nVertLevels + block_ptr->nVertLevels; i++)
						fprintf(stderr, "%f ", block_ptr->uVertexVelocities[1][i]);

				}
			
				// dprint("pred particles %ld, rank %d", block_ptr->particles.size(), block_ptr->gid);

				// dprint("gid %d particles %ld", world.rank(), block_ptr->particles.size());
				// send sampled seeds to origin processes
				send_sampled_seeds_to_origin(master_pred, assigner);

				// dprint("gid %d particles %ld", world.rank(), block_ptr->particles.size());

				
				// streamline advection
				dbgmsg = "###particles before pred";
				print_global_num_particles(dbgmsg, world, master_pred, framenum);

				

				master_pred.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {
					streamline sl(*b, dtSim, dtParticle);

					

					bool val = trace_particles(b,
												icp,
												assigner,
												max_steps,
												&sl,
												true,
												nsteps,
												particles_hold,
												skip_rate);

					return val;
				});

				dbgmsg = "###particles after pred";
				print_global_num_particles(dbgmsg, world, master_pred, framenum);


				// save block and initial links
				// block_ptr = reinterpret_cast<block*> (master_pred.release(master_pred.lid(gids[0])));
				
				
				// move particles_hold to master_baladv->particles
				master_baladv.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					b->particles = std::move(particles_hold);
				});		

				dbgmsg = "###particles after hold brought back";
				print_global_num_particles(dbgmsg, world, master_baladv, framenum);

				// kd-tree based partition
				// load based repartition
				bool wrap = false;
				size_t samples = 512;
				diy::kdtree(master_baladv, assigner, DIM, domain, &block::particles, samples, wrap);


				dbgmsg = "### particles before pos data and seeds after bal";
				print_global_num_particles(dbgmsg, world, master_baladv, framenum);

				// process data and seeds, also update local_gcIds, and gcIdToGid
				master_baladv.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					b->position_data_and_seeds_after_balance(world, particles_hold);
				});	

				dbgmsg = "### particles after pos data and seeds after bal";
				print_global_num_particles(dbgmsg, world, master_baladv, framenum);


			}else{

				// sort data based on current partition after moving latest to regular position
				sort_data(master_baladv, assigner);

			}

			

			
			if (framenum > 1){
				if (world.rank()==0)
					dprint("starting pathline advection frame %d", framenum);

				ghost_exchange(master_baladv, assigner, false);
				

				// interpolate to vertex
				interpolate_cell_to_vertex(master_baladv, false);

				dbgmsg = "###particles before adv";
				print_global_num_particles(dbgmsg, world, master_baladv, framenum);	

				// // pathline advection
				// master_baladv.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {
				// 	pathline pl(*b, dtSim, dtParticle);

				// 	bool val = trace_particles(b,
				// 								icp,
				// 								assigner,
				// 								max_steps,
				// 								&pl,
				// 								true,
				// 								nsteps,
				// 								particles_hold,
				// 								skip_rate);

				// 	return val;
				// });

				dbgmsg = "###particles after adv";
				print_global_num_particles(dbgmsg, world, master_baladv, framenum);	


				
			}
			
			
			particles_hold.clear();

		}
			
			
		

	}		




	if (world.rank()==0){
		dprint("done ");
	}
	
}

int main(int argc, char* argv[])
{
	
	// dprint("In con");

	block mpas1;
	double dtSim = 7200, dtParticle = 300;
	pathline pl(mpas1, dtSim, dtParticle);

	// int nblocks = 0;
	// int threads = 2;
	// int mem_blocks = -1;
	int max_steps = 5;
	int check = 0;                  // write out traces
	int gid;

	int ndims = 3;           // domain dimensions
	string particle_file;    // input file name
	string gp_file;    		// input file name
	// double dtSim = 5000*7200*10;
	// double dtParticle = 300000;
	size_t nCells;
	int pred_percent = 10;
	int seed_rate = 100;
	int skip_rate = 10;
	bool prediction = true;
	int rebal_interval = 1;


	// define the workflow
	Workflow workflow;
	//make_wflow(workflow);
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");

	MPI_Init(NULL, NULL);

	// create decaf
	Decaf *decaf = new Decaf(MPI_COMM_WORLD, workflow);

	diy::mpi::communicator world(decaf->con_comm_handle());
	// dprint("diy world size: %d", world.size() );

	// if (world.rank() == 1)
	//    {

	//        volatile  int dwait=0;
	//        fprintf(stderr , "pid %ld  waiting  for  debugger\n"
	//            , (long)getpid ());
	//            while(dwait==0) { /*  change  ’i’ in the  debugger  */ }
	//    }
	//    world.barrier();

	diy::FileStorage storage("./DIY.XXXXXX");
	int nblocks = world.size();
	int threads = 1;
	int mem_blocks = -1;

	using namespace opts;

	 // command-line ags
    Options ops(argc, argv);

	ops >> Option('b', "blocks", nblocks, "Total number of blocks to use")
		>> Option('c', "check", check, "Write out traces for checking");

	 if (ops >> Present('h', "help", "show help") ||
        !(ops >> PosOption(particle_file) >> PosOption(gp_file) >> PosOption(prediction) >> PosOption(dtSim) >> PosOption(dtParticle)
				>> PosOption(seed_rate) >> PosOption(pred_percent)>> PosOption(skip_rate) >> PosOption(rebal_interval)))
    {
        if (world.rank() == 0)
        {
            fprintf(stderr, "Check ops usage \n", argv[0]);
            cout << ops;
        }
        return 1;
    }


	// diy::Master master(world,
	// 				   threads,
	// 				   mem_blocks,
	// 				   NULL,
	// 				   NULL,
	// 				   &storage,
	// 				   NULL,
	// 				   NULL);

	// diy::RoundRobinAssigner assigner(world.size(), nblocks);

	// std::vector<int> gids;					 // global ids of local blocks
	// assigner.local_gids(world.rank(), gids); // get the gids of local blocks
	// mpas1.gid = gids[0];

	// diy::RoundRobinAssigner pred_assigner(world.size(), nblocks);
	// pred_assigner.local_gids(world.rank(), gids); // get the gids of local blocks
	// mpas1.gid = gids[0];


	// dprint("max_steps %d skipt rate %d", max_steps, skip_rate);
	con(decaf, mpas1, pl, gp_file, dtSim, dtParticle, particle_file, seed_rate, pred_percent, max_steps, skip_rate, prediction, nblocks, rebal_interval);

	delete decaf;

	// write out stats
	size_t nsteps_global=0;
	diy::mpi::reduce(world, nsteps, nsteps_global, 0, std::plus<size_t>());

	size_t maxsteps_global=0;
    diy::mpi::reduce(world, nsteps, maxsteps_global, 0, diy::mpi::maximum<size_t>());

	if (world.rank()==0)
		fprintf(stderr, "predd , %d, nsteps_global , %ld, maxsteps_global , %ld, time_predrun , %f, time_final, %f, time_readdata, %f, worldsize, %d, minsteps, %ld, dtSim %f, dtParticle, %f,\n", prediction, nsteps_global, maxsteps_global, time_predrun, time_final, time_readdata, world.size(), 0, dtSim, dtParticle);


	// MPI_Finalize(); // called by diy for consumers
}
