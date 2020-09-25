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

void enq_halo_info(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner)
{

	Halo h;

	std::vector<Halo> halos = b->get_halo_info();
	for (size_t i = 0; i < halos.size(); i++)
	{

		int dest_gid = halos[i].dest_gid;
		int dest_proc = assigner.rank(dest_gid);
		diy::BlockID dest_block = {dest_gid, dest_proc};
		cp.enqueue(dest_block, halos[i]);
	}
}

void process_halo_req(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int framenum)
{

	vector<int> in;
	cp.incoming(in);
	for (int i = 0; i < in.size(); i++)
	{
		if (cp.incoming(in[i]).buffer.size() > 0)
		{
			Halo incoming_halo;
			cp.dequeue(in[i], incoming_halo);
			b->process_halo_req(incoming_halo, framenum);

			int dest_gid = incoming_halo.src_gid;
			int dest_proc = assigner.rank(dest_gid);
			diy::BlockID dest_block = {dest_gid, dest_proc};
			cp.enqueue(dest_block, incoming_halo);
		}
	}
}

void deq_halo_info(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int framenum)
{

	vector<int> in;
	cp.incoming(in);
	for (int i = 0; i < in.size(); i++)
	{
		if (cp.incoming(in[i]).buffer.size() > 0)
		{
			Halo incoming_halo;
			cp.dequeue(in[i], incoming_halo);
			b->update_halo_info(incoming_halo, framenum);
		}
	}
}

// for dynamic processing and queuing
void process_halo_dynamic(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int framenum, block)
{

	b->process_halo_dynamic(framenum);

	// for each destination in b->halo_info
	for (size_t i = 0; i < b->halo_info.size(); i++)
	{

		int dest_gid = b->halo_info[i].src_gid;
		int dest_proc = assigner.rank(dest_gid);
		diy::BlockID dest_block = {dest_gid, dest_proc};
		cp.enqueue(dest_block, b->halo_info[i]);
	}
}

void deq_halo_dynamic(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int framenum)
{

	vector<int> in;
	cp.incoming(in);
	for (int i = 0; i < in.size(); i++)
	{
		if (cp.incoming(in[i]).buffer.size() > 0)
		{
			Halo incoming_halo;
			cp.dequeue(in[i], incoming_halo);
			b->update_halo_dynamic(incoming_halo, framenum);
		}
	}
}

void all_gather_data(diy::mpi::communicator &world, std::vector<int> &data_bar_int, std::vector<double> &data_bar_dbl, int data_id)
{

	if (data_bar_int.size() > 1)
	{
		std::vector<std::vector<int>> out;
		diy::mpi::all_gather(world, data_bar_int, out);
		data_bar_int.clear();

		// if (data_id==6){
		// 	int cur_start = 0;
		// 	for (size_t i=0; i<out.size(); i++){

		// 		for (int &d: out[i]){
		// 				// dprint("%d ",d);
		// 			d += cur_start;
		// 		}
		// 		// exit(0);

		// 		cur_start += out[i].size();
		// 		dprint("cur_start %d", cur_start);
		// 	}
		// }

		for (size_t i = 0; i < out.size(); i++)
		{
			data_bar_int.insert(data_bar_int.end(), out[i].begin(), out[i].end());
		}
		// if (data_id==6){
		// 	for (size_t i =0 ; i<data_bar_int.size(); i++)
		// 		if (data_bar_int[i] == 376)
		// 			dprint("376 found at %ld", i);
		// }
	}

	if (data_bar_dbl.size() > 1)
	{
		std::vector<std::vector<double>> out;
		diy::mpi::all_gather(world, data_bar_dbl, out);
		data_bar_dbl.clear();

		for (size_t i = 0; i < out.size(); i++)
		{
			data_bar_dbl.insert(data_bar_dbl.end(), out[i].begin(), out[i].end());
		}
	}
}

bool first_done = false;

// consumer
void con(Decaf *decaf, diy::Master &master, diy::RoundRobinAssigner &assigner, block &mpas1, pathline &pl, std::string &fname_decoupled_particles)
{

	// mpaso mpas_c;
	vector<pConstructData> in_data;
	diy::mpi::communicator world(decaf->con_comm_handle());

	// recomputed particle file ("particles2.nc") includes global cell indices for seeds
	// std::string fname_particles = "particles2.nc", fname_graphinfo = "graph.info.part." + std::to_string(world.size());
	std::string fname_particles = fname_decoupled_particles, fname_graphinfo = "graph.info.part." + std::to_string(world.size());

	// // build kd tree
	// mpas1.create_cells_kdtree();
	dprint("In con");

	

	double time_trace = 0, time_gatherd = 0, time_lb = 0;

	while (decaf->get(in_data))
	{	

		
		int n_velocityX = 0;
		std::vector<int> indexToCellID;
		std::vector<double> data_bar_dbl;
		std::vector<int> data_bar_int;
		int fir_cell_idx = 999;

		Bounds domain{DIM};
		for (unsigned i = 0; i < DIM; ++i)
		{
			domain.min[i] = -6371329.; // 6371229.
			domain.max[i] = 6371329.;
		}

		// get the values and add them
		for (size_t i = 0; i < in_data.size(); i++)
		{

			SimpleFieldi d_metadata = in_data[i]->getFieldData<SimpleFieldi>("field_id");
			SimpleFieldi frame_num = in_data[i]->getFieldData<SimpleFieldi>("frame_num");
			data_bar_int = in_data[i]->getFieldData<VectorFieldi>("data_i").getVector();
			data_bar_dbl = in_data[i]->getFieldData<VectorFliedd>("data_d").getVector();

			// printf("fieldid %d \n", d_metadata.getData());

			int data_id = d_metadata.getData();
			int framenum = frame_num.getData();

			// dprint("data_id %d", data_id);

			

			if (framenum == 1 && first_done == true)
			{
				// to deal with the output stream for restart files
				// dprint("BREAKING");
				break;
			}
			if (framenum > 1)
			{
				first_done = true;
			}

			// all_gather data_bar
			world.barrier();
			double time_0 = MPI_Wtime();

			all_gather_data(world, data_bar_int, data_bar_dbl, data_id);
			
			world.barrier();
			double time_1 = MPI_Wtime();
			time_gatherd += time_1 - time_0;

			mpas1.update_data(data_id, framenum, data_bar_int, data_bar_dbl);

			// if (data_id==17){
			// 	for (size_t i=0; i<mpas1.vertVelocityTop[framenum%2].size();i++)
			// 		if (std::abs(mpas1.vertVelocityTop[framenum%2][i]-0.000000372406)<0.00000000001)
			// 			dprint("found at %ld, %.10f", i, mpas1.vertVelocityTop[framenum%2][i]);

			// }

			if (0)
			if (framenum == 1 && data_id == 15)
			{ // all static arrived

				// create links
				// diy::Link* link = new diy::Link;
				RCLink *l = new RCLink(DIM, domain, domain);
				// add links and block to master
				master.add(mpas1.gid, &mpas1, l);

				// initialize particles
				//if (world.rank() == 0)
				//{
					// mpas1.generate_new_particle_file();
					// mpas1.init_seeds_mpas(fname_particles, framenum, world.rank());
				//}
			}
			else if (framenum > 1 && data_id == 13)
			{   // all dynamic arrived
				// if ((data_id==13 && framenum>1) || (data_id==15 && framenum==1)){ // all dynamic arrived

				// kdtre based balance

				dprint("IN HERE %d, mpas1.velXv %d, %f", framenum, mpas1.velocityXv[0].size(), mpas1.velocityXv[0][0]); //mpas1.velocityXv[0][0]

				bool wrap = false;
				size_t samples = 512;

				world.barrier();
				double time_0 = MPI_Wtime();

				// ### Uncomment following line to enable kdtree ###
				// diy::kdtree(master, assigner, DIM, domain, &block::particles, samples, wrap);

				// master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
				// {
				// 	dprint("inside");
				// 	dprint ("psize %ld", b->particles.size());
				// });

				// advect

				world.barrier();
				double time_1 = MPI_Wtime();
				time_lb += time_1 - time_0;

				//do{
				master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					// update field
					pl.update_velocity_vectors(*b, framenum);

					// dprint("here2");
					// trace particles
					// pl.compute_epoch(b, framenum);
					// dprint("here3 rank %d", world.rank());

					// b->parallel_write_simstep_segments(world, framenum);
				});

				world.barrier();
				double time_2 = MPI_Wtime();
				
				time_trace += time_2 - time_1;
				//}while(reduce of all remaining particles not zero)
				if (world.rank() == 0)
					dprint("Completed dynamic frame %d", framenum);
				
			}

		
	
		} // looping through all in_data items

		


	} // outer while loop decaf get

	if (world.rank()==0)
		dprint("rank, %d, times trace lb gatherd, %f, %f, %f, ", world.size(), time_trace, time_lb, time_gatherd);

	// dprint("writing segments");
	// write segments to file in parallel
	// master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
	// dprint("writing");
	// b->parallel_write_segments(world, 0);
	// });
}

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
					 std::vector<EndPt>& particles_hold, 
					 int skip_rate){
				bool val = true;

				// b->particles_store.clear();
				deq_incoming_iexchange(b, cp);
				// dprint("calling trace_particles");
				val = pl.compute_streamlines(b, cp, assigner, prediction, nsteps, particles_hold, skip_rate);

				return val;

}

void print_block(block* b, const diy::Master::ProxyWithLink& cp, bool verbose, int worldsize)
{
  RCLink*  link      = static_cast<RCLink*>(cp.link());

  fmt::print("bounds,{},: [,{},{},{},] - [,{},{},{},] ({} neighbors): {} points, world, {},\n",
                  cp.gid(),
                  link->bounds().min[0], link->bounds().min[1], link->bounds().min[2],
                  link->bounds().max[0], link->bounds().max[1], link->bounds().max[2],
                  link->size(), b->particles.size(), worldsize);

//   for (int i = 0; i < link->size(); ++i)
//   {
//       fmt::print("  ({},{},({},{},{})):",
//                       link->target(i).gid, link->target(i).proc,
//                       link->direction(i)[0],
//                       link->direction(i)[1],
//                       link->direction(i)[2]);
//       const Bounds& bounds = link->bounds(i);
//       fmt::print(" [{},{},{}] - [{},{},{}]\n",
//               bounds.min[0], bounds.min[1], bounds.min[2],
//               bounds.max[0], bounds.max[1], bounds.max[2]);
//   }

//   if (verbose)
//     for (size_t i = 0; i < b->points.size(); ++i)
//       fmt::print("  {} {} {}\n", b->points[i][0], b->points[i][1], b->points[i][2]);
}


int main(int argc, char* argv[])
{
	dprint("starting main....");

	MPI_Init(NULL, NULL);

	// define the workflow
	Workflow workflow;
	//make_wflow(workflow);
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");
	// create decaf
	Decaf *decaf = new Decaf(MPI_COMM_WORLD, workflow);
	diy::mpi::communicator world(decaf->con_comm_handle());


	diy::FileStorage storage("./DIY.XXXXXX");
	int nblocks = world.size();
	int threads = 2;
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
	int pred_percent = 10;
	int seed_rate = 100;
	int skip_rate = 10;

	bool prediction = true;
	double time_prep=0, time_predrun=0, time_kdtree=0, time_readdata=0, time_filter=0, time_final=0, time_predrun_loc=0, time_final_loc, time_trace = 0;
	std::atomic<bool> done{false};
	std::vector<size_t> steps_per_interval;
	std::vector<size_t> particles_in_core;
	size_t np_core = 0;
	std::vector<int> gcIdxToGid;



	block mpas1;
	// double dtSim = 7200, dtParticle = 300;
	pathline pl(mpas1, dtSim, dtParticle);

	
	

	
	// dprint("diy world size: %d", world.size() );

	// if (world.rank() == 1)
	//    {

	//        volatile  int dwait=0;
	//        fprintf(stderr , "pid %ld  waiting  for  debugger\n"
	//            , (long)getpid ());
	//            while(dwait==0) { /*  change  ’i’ in the  debugger  */ }
	//    }
	//    world.barrier();

	// diy::FileStorage storage("./DIY.XXXXXX");
	// int nblocks = world.size();
	// int threads = 1;
	// int mem_blocks = -1;



	using namespace opts;

	 // command-line ags
    Options ops(argc, argv);

	ops >> Option('b', "blocks", nblocks, "Total number of blocks to use")
		>> Option('c', "check", check, "Write out traces for checking");

	 if (ops >> Present('h', "help", "show help") ||
        !(ops >> PosOption(particle_file) >> PosOption(gp_file) >> PosOption(prediction) >> PosOption(dtSim) >> PosOption(dtParticle)
				>> PosOption(seed_rate) >> PosOption(pred_percent)>> PosOption(skip_rate)))
    {
        if (world.rank() == 0)
        {
            fprintf(stderr, "Check ops usage \n", argv[0]);
            cout << ops;
        }
        return 1;
    }

	dprint("%s", particle_file.c_str());
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

	std::vector<int> gids;					 // global ids of local blocks
	assigner.local_gids(world.rank(), gids); // get the gids of local blocks
	mpas1.gid = gids[0];


	
	con(decaf, master, assigner, mpas1, pl, particle_file);


	// start pred


	// // add master to block 
	// std::vector<int> gids;                     // global ids of local blocks
    // assigner.local_gids(world.rank(), gids);   // get the gids of local blocks

	for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
    {	block* b = new block();
		int gid = gids[i];
		b->gid = gid;

		diy::Link*    link = new diy::Link;  // link is this block's neighborhood


		// read mpas data
		std::string fname_data = "output/output.nc";
		b->loadMeshFromNetCDF_CANGA(world, fname_data, 0);
		

		
		// read and add link
		std::string fname_graph = "graph.info", fname_graphpart = gp_file;
		//"graph.info.part." + std::to_string(world.size());
		dprint("fname_graph %s, fname_graphpart %s", fname_graph.c_str(), fname_graphpart.c_str());

        set<int> links;
        b->create_links(fname_graph, fname_graphpart, links);
		dprint("rank %d links %ld",world.rank(), links.size());
		for (int lgid: links){
			dprint("rank %d links to %ld", world.rank(), i);
			diy::BlockID  neighbor;
			neighbor.gid  = lgid;                     // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			link->add_neighbor(neighbor);
		}


	
		// init particles
		std::string fname_particles = particle_file; //"particles.nc";
		b->init_seeds_particles(world, fname_particles, 0, seed_rate);

		

		b->init_partitions(); // must be called after loading mpas data, after setting b->gid, and after init_seeds_partitions
		nCells = b->gcIdxToGid.size();
		

		master.add(gid, b, link);
	}

	

	std::vector<EndPt> particles_hold; // staging area for holding to-be-predicted particles
    if (0)
	if (prediction)
	{	
		Bounds domain{DIM};
		for (unsigned i = 0; i < DIM; ++i)
		{
			domain.min[i] = -6371329. - 100; // 6371229.
			domain.max[i] = 6371329. + 100;
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
			std::mt19937 g(6);
			std::shuffle(b->particles.begin(), b->particles.end(), g);

			// size_t pred_size = b->particles.size()/10;
			size_t pred_size = static_cast<size_t>(floor((float(pred_percent)/100)*b->particles.size()));

			// dprint("pred_size %d", pred_size);

			// b->particles_hold.insert(std::end(b->particles_hold), std::begin(b->particles) + pred_size, std::end(b->particles));
			// particles_hold.insert(std::end(particles_hold), std::begin(b->particles) + 0, std::end(b->particles));
			particles_hold.insert(std::end(particles_hold), std::begin(b->particles) + pred_size, std::end(b->particles));
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

		

		// prediction advection using iexchange

		master.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {
			return true;
		});
		
		master.iexchange([&](block *b, const diy::Master::ProxyWithLink &icp) -> bool {

			pathline pl(*b, dtSim, dtParticle);

			// update_velocity_vectors: both timesteps point to same vectors
			pl.set_velocity_vectors(*b);

			// dprint("before trace");

			// bool val = trace_particles(b,
			// 						   icp,
			// 						   assigner,
			// 						   max_steps,
			// 						   pl,
			// 						   true,
			// 						   nsteps,
			// 						   particles_hold, 
			// 						   skip_rate);
			bool val = true;
			dprint("after tracee %d", val);

			return val;
		});
		dprint("before third barrier");
		world.barrier();
		double time2 = MPI_Wtime(); 
		time_predrun = time2 - time1;

		dprint("after third barrier");


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
				// if (cell_cen.glCellIdx==196056 || cell_cen.glCellIdx==48616){
				// 	dprint("cell_cen.glCellIdx %d starts from gid %d", cell_cen.glCellIdx , b->gid);
				// }
			}

			// b->particles = std::move(b->particles_hold);
			// dprint("particles size %ld", b->particles.size());
			block_ptr = std::move(b);
		});

		master_kdt.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			b->particles = std::move(particles_hold);
		});


		world.barrier();
		double time3 = MPI_Wtime();
		time_prep += time3 - time3;
		
		// load based repartition
		bool wrap = false;
		size_t samples = 512;
		diy::kdtree(master_kdt, assigner, DIM, domain, &block::particles, samples, wrap);

		bool verbose = false;
		master_kdt.foreach([&](block* b, const diy::Master::ProxyWithLink& cp) { print_block(b,cp,verbose, world.size()); });

		master_kdt.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
			particles_hold = std::move(b->particles);
			gid = b->gid;
		});

		world.barrier();
		double time4 = MPI_Wtime();
		time_kdtree = time4 - time3;

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
		// nsteps = 0;

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

			// if (b->gid == 1138)
			// {
			// 	diy::BlockID neighbor;
			// 	neighbor.gid = 1150;						 // gid of the neighbor block
			// 	neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			// 	link->add_neighbor(neighbor);
			// }
			// if (b->gid == 1599)
			// {
			// 	diy::BlockID neighbor;
			// 	neighbor.gid = 1623;						 // gid of the neighbor block
			// 	neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			// 	link->add_neighbor(neighbor);
			// }

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
					// if (b->particles[i].glCellIdx==196056 || b->particles[i].glCellIdx==48616 ){
					// 	dprint("cell_cen.glCellIdx %d goes to gid %d", b->particles[i].glCellIdx , b->gid);
					// }
				}else if (b->particles[i].predonly==0){
					b->particles[i].pt = b->particles[i].pt_hold;
					particles_hold.push_back(b->particles[i]);
				}
			}

			

			// dprint("post partition %d, particles %ld", postsum, b->particles.size());
			b->particles = std::move(particles_hold);
		});

		dprint("Starting second advection .............. ");


			world.barrier();
        	double time6 = MPI_Wtime();	
			time_filter = time6 - time4;

                

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
											particles_hold,
											skip_rate);

				

				
				return val;
			});

			world.barrier();
			double time7 = MPI_Wtime();
			time_final = time7 - time6;

			dprint("done final advection");

			master_final.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {	
				gcIdxToGid = std::move(b->gcIdxToGid);
			});

	}



	// end pred


	delete decaf;

	// MPI_Finalize(); // called by diy for consumers
}
