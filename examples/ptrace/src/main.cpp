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
void con(Decaf *decaf, diy::Master &master, diy::RoundRobinAssigner &assigner, block &mpas1, pathline &pl)
{

	// mpaso mpas_c;
	vector<pConstructData> in_data;
	diy::mpi::communicator world(decaf->con_comm_handle());

	// recomputed particle file ("particles2.nc") includes global cell indices for seeds
	std::string fname_particles = "particles2.nc", fname_graphinfo = "graph.info.part." + std::to_string(world.size());

	// build kd tree
	mpas1.create_cells_kdtree();

	

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

			if (framenum == 1 && data_id == 15)
			{ // all static arrived

				// create links
				// diy::Link* link = new diy::Link;
				RCLink *l = new RCLink(DIM, domain, domain);
				// add links and block to master
				master.add(mpas1.gid, &mpas1, l);

				// initialize particles
				if (world.rank() == 0)
				{
					mpas1.generate_new_particle_file();
					mpas1.init_seeds_mpas(fname_particles, framenum, world.rank());
				}
			}
			else if (framenum > 1 && data_id == 13)
			{   // all dynamic arrived
				// if ((data_id==13 && framenum>1) || (data_id==15 && framenum==1)){ // all dynamic arrived

				// kdtre based balance

				// dprint("IN HERE %d", framenum);

				bool wrap = false;
				size_t samples = 512;

				world.barrier();
				double time_0 = MPI_Wtime();

				diy::kdtree(master, assigner, DIM, domain, &block::particles, samples, wrap);
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
					pl.compute_epoch(b, framenum);
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

int main()
{

	block mpas1;
	double dtSim = 7200, dtParticle = 300;
	pathline pl(mpas1, dtSim, dtParticle);

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

	con(decaf, master, assigner, mpas1, pl);

	delete decaf;

	// MPI_Finalize(); // called by diy for consumers
}
