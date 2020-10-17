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

typedef diy::RegularContinuousLink RCLink;
typedef diy::ContinuousBounds Bounds;
static const unsigned DIM = 3;

void copy_cell_data_to_endpt(EndPt &cell_cen, block &b, int idx){

	cell_cen.cell_data.resize(6);	
	cell_cen.cell_data[0].clear();
	cell_cen.cell_data[0].insert(cell_cen.cell_data[0].end(), b.velocityX.begin()+idx*b.nVertLevels, b.velocityX.begin()+idx*b.nVertLevels + b.nVertLevels);
	cell_cen.cell_data[1].clear();
	cell_cen.cell_data[1].insert(cell_cen.cell_data[1].end(), b.velocityY.begin()+idx*b.nVertLevels, b.velocityY.begin()+idx*b.nVertLevels + b.nVertLevels);
	cell_cen.cell_data[2].clear();
	cell_cen.cell_data[2].insert(cell_cen.cell_data[2].end(), b.velocityZ.begin()+idx*b.nVertLevels, b.velocityZ.begin()+idx*b.nVertLevels + b.nVertLevels);
	cell_cen.cell_data[3].clear();
	cell_cen.cell_data[3].insert(cell_cen.cell_data[3].end(), b.vertVelocityTop.begin()+idx*b.nVertLevels, b.vertVelocityTop.begin()+idx*b.nVertLevels + b.nVertLevels);
	cell_cen.cell_data[4].clear();
	cell_cen.cell_data[4].insert(cell_cen.cell_data[4].end(), b.zTop.begin()+idx*b.nVertLevels, b.zTop.begin()+idx*b.nVertLevels + b.nVertLevels);
	cell_cen.cell_data[5].clear();
	cell_cen.cell_data[5].insert(cell_cen.cell_data[5].end(), b.zMid.begin()+idx*b.nVertLevels, b.zMid.begin()+idx*b.nVertLevels + b.nVertLevels);

}

void copy_endpt_to_cell_data(EndPt &cell_cen, block &b){

	int idx = cell_cen.glCellIdx;
	// b.velocityX.clear();	
	b.velocityX.insert(b.velocityX.begin()+idx*b.nVertLevels, cell_cen.cell_data[0].begin(), cell_cen.cell_data[0].end());
	// b.velocityY.clear();
	b.velocityY.insert(b.velocityY.begin()+idx*b.nVertLevels, cell_cen.cell_data[1].begin(), cell_cen.cell_data[1].end());
	// b.velocityZ.clear();
	b.velocityZ.insert(b.velocityZ.begin()+idx*b.nVertLevels, cell_cen.cell_data[2].begin(), cell_cen.cell_data[2].end());
	// b.vertVelocityTop.clear();
	b.vertVelocityTop.insert(b.vertVelocityTop.begin()+idx*b.nVertLevels, cell_cen.cell_data[3].begin(), cell_cen.cell_data[3].end());
	// b.zTop.clear();
	b.zTop.insert(b.zTop.begin()+idx*b.nVertLevels, cell_cen.cell_data[4].begin(), cell_cen.cell_data[4].end());
	// b.zMid.clear();
	b.zMid.insert(b.zMid.begin()+idx*b.nVertLevels, cell_cen.cell_data[5].begin(), cell_cen.cell_data[5].end());
	
}


void enqueue_ghost_cells_to_nbrs(std::vector<int> &local_gcIds, block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, std::vector<int> &gcIdxToGid_global, block &mpas1){

	
	
	std::map<diy::BlockID, std::vector<EndPt>>   outgoing_cells;
	for (size_t i = 1; i < local_gcIds.size(); i++)
	{	
		int gcId = local_gcIds[i];
		// if (i<10)
		// 	dprint("i %d size %ld", i, b->cell_nbrs[gcId].size())
		for (size_t j=0; j<b->cell_nbrs[gcId].size();j++){
			int nbr_gcId = b->cell_nbrs[gcId][j]; 
			// dprint("nbr_gcId %d, %d %ld", nbr_gcId, gcId, j);
			int nbr_gid = b->gcIdToGid[b->cell_nbrs[gcId][j]];
			if (nbr_gid != mpas1.gid){
				// send info of gcId to nbr_gid, store in cell_data

				EndPt endpt;
				endpt.glCellIdx = gcId - 1;
				copy_cell_data_to_endpt(endpt, *b, endpt.glCellIdx);

				int dest_gid = nbr_gid;
				int dest_proc = assigner.rank(dest_gid);
				diy::BlockID dest_block = {dest_gid, dest_proc};

				outgoing_cells[dest_block].push_back(endpt);
				// if (b->gid==2){
				// 	dprint("gcId %d, nbr %d, nbr_gcId %d nbr_gid %d", gcId, dest_gid, nbr_gcId, nbr_gid);
				// }

			}
		}

	
		
	}

	for (auto elem: outgoing_cells){
		cp.enqueue(elem.first, elem.second);
	}
	
	dprint("gid %d outsize %ld", b->gid, outgoing_cells.size());
	// if (b->gid==2){
	// 	for (auto i : outgoing_cells){
	// 		dprint("outs %d", i.second.size());
	// 	}

	// }
}

void deque_ghost_cells_from_nbrs(block *b, const diy::Master::ProxyWithLink &cp){

	vector<int> in;
	cp.incoming(in);
	for (int i = 0; i < in.size(); i++)
	{
		if (cp.incoming(in[i]).buffer.size() > 0)
		{
			std::vector<EndPt> incoming_endpts;
			cp.dequeue(in[i], incoming_endpts);
			// b->process_halo_req(incoming_halo, framenum);
			dprint("cp.gid %d, numdeq %d", cp.gid(), incoming_endpts.size());

			for (size_t j=0; j<incoming_endpts.size(); j++){
				copy_endpt_to_cell_data(incoming_endpts[j], *b);
			}
		}
	}
	
	
}

bool first_done = false;

// consumer
void con(Decaf *decaf, diy::Master &master, diy::RoundRobinAssigner &assigner, block &mpas1, pathline &pl, string &gp_file, double dtSim, double dtParticle)
{

	// mpaso mpas_c;
	vector<pConstructData> in_data;
	diy::mpi::communicator world(decaf->con_comm_handle());

	// recomputed particle file ("particles2.nc") includes global cell indices for seeds
	std::string fname_particles = "particles2.nc", fname_graphinfo = "graph.info.part." + std::to_string(world.size());

	// build kd tree
	// mpas1.create_cells_kdtree();

	// create block, read mpas data and add to master
	// add master to block 
	std::vector<int> gids;                     // global ids of local blocks
    assigner.local_gids(world.rank(), gids);   // get the gids of local blocks

	for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
    {	block* b = new block();
		int gid = gids[i];
		b->gid = gid;

		diy::Link*    link = new diy::Link;  // link is this block's neighborhood

		std::string fname_graph = "graph.info", fname_graphpart = gp_file;
		//"graph.info.part." + std::to_string(world.size());
		dprint("fname_graphpart %s", fname_graphpart.c_str());
		// read mpas data
		std::string fname_data = "output2.nc";
		b->loadMeshFromNetCDF_CANGA(world, fname_data, 0, fname_graph);
		

		// read and add link
        set<int> links;
        b->create_links(fname_graph, fname_graphpart, links);
		// dprint("rank %d links %ld", gid, links.size());
		if(gid==2){
			for (auto i:links)
				dprint("link %d", i);
		}
		
			
		for (int lgid: links){
			// dprint("rank %d links to %ld", world.rank(), i);
			diy::BlockID  neighbor;
			neighbor.gid  = lgid;                     // gid of the neighbor block
            neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			link->add_neighbor(neighbor);
		}
	

		master.add(gid, b, link);
	}


	// do initial ghost exchange





	double time_trace = 0, time_gatherd = 0, time_lb = 0;

	// int cnt = 0;

	while (decaf->get(in_data))
	{	
		
		// continue;

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
			// SimpleFieldi d_metadata = in_data.at("in")->getFieldData<SimpleFieldi>("field_id");
			SimpleFieldi frame_num = in_data[i]->getFieldData<SimpleFieldi>("frame_num");
			data_bar_int = in_data[i]->getFieldData<VectorFieldi>("data_i").getVector();
			data_bar_dbl = in_data[i]->getFieldData<VectorFliedd>("data_d").getVector();

			// printf("fieldid %d \n", d_metadata.getData());
			
				int data_id = d_metadata.getData();
				int framenum = frame_num.getData();
				if (world.rank()==0){
					dprint("data_id %d, framenum %ld, %ld %ld, %d %f", data_id, framenum, data_bar_int.size(), data_bar_dbl.size(), data_bar_int[0], data_bar_dbl[0]);
				}
			
			mpas1.update_data(data_id, framenum, data_bar_int, data_bar_dbl);


		

		
		
	
		

		if (data_id == 15){ // all static data has arrived

			std::vector<EndPt> particles_hold; // staging area for holding to-be-predicted particles
			Bounds domain{DIM};
			for (unsigned i = 0; i < DIM; ++i)
			{
				domain.min[i] = -6371329. - 100; // 6371229. for EC
				domain.max[i] = 6371329. + 100;
			}
				

			diy::FileStorage storage("./DIY.XXXXXX");
			int threads = 2;
			int mem_blocks = -1;
			diy::Master master_kdt(world,
						threads,
						mem_blocks,
						NULL,
						NULL,
						&storage,
						NULL,
						NULL);

			for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
				{	block* b = new block();
					int gid = gids[i];
					b->gid = gid;

					RCLink*         l   = new RCLink(DIM, domain, domain);
					master_kdt.add(gid, b, l);

				}

		


				// select prediction particles; store tag along particles in b->particles_hold
				master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					std::random_device rd;
					std::mt19937 g(6);
					std::shuffle(b->particles.begin(), b->particles.end(), g);

					// size_t pred_size = static_cast<size_t>(floor((float(pred_percent)/100)*b->particles.size()));

					// particles_hold.insert(std::end(particles_hold), std::begin(b->particles) + pred_size, std::end(b->particles));
					// b->particles.resize(pred_size);

					for (EndPt& ep: particles_hold){ // updating p{xyz} to cell{xyz} for load based sorting
						ep.pt_hold = ep.pt;
						ep[0] = b->xCell[ep.glCellIdx];
						ep[1] = b->yCell[ep.glCellIdx];
						ep[2] = b->zCell[ep.glCellIdx];
					}

				});

				// preparing to compute links by computing gcIdxToGid_global	
				std::vector<int> gcIdToGid_local(mpas1.nCells+1);
				std::vector<int> gcIdToGid_global(mpas1.nCells+1);
				std::vector<int> local_gcIds; // stores the list of gcIdx in currect process

				// replace particle p{xyz} with cell{xyz}, move b->particles_hold to kdt master's b->particles
				master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {

					dprint("mpas1.indexToCellID %ld, mpas1.velocityX %ld", mpas1.indexToCellID.size(), mpas1.velocityX.size());

					int idx = 0;
					for (auto i: mpas1.indexToCellID){
						EndPt cell_cen;
						cell_cen[0] = b->xCell[i];
						cell_cen[1] = b->yCell[i];
						cell_cen[2] = b->zCell[i];
						cell_cen.predonly = 2; // indicates this is the cell center point
						cell_cen.glCellIdx = i - 1;
						cell_cen.source_gid = b->gid;
						local_gcIds.push_back(i);
						gcIdToGid_local[i] = b->gid; // needs to be filled before reduction and ghost exchange

						// copy cell velocities from mpas1 to cell_cen (mpas1 -> cell_cen)
						copy_cell_data_to_endpt(cell_cen, mpas1, idx);
						
						// copy cell_cen to b for prediction advection
						copy_endpt_to_cell_data(cell_cen, *b);

						// push cell_cen data to particles_hold for partitioning after prediction
						particles_hold.push_back(cell_cen);

						idx++;

						

					}


				});

				diy::mpi::all_reduce (world, gcIdToGid_local, gcIdToGid_global,  std::plus<int>());

				// do initial ghost exchange
				master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {	
					enqueue_ghost_cells_to_nbrs(local_gcIds, b, cp, assigner, gcIdToGid_global, mpas1);
				});
				master.exchange();
				master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					deque_ghost_cells_from_nbrs(b, cp);
					
					// compute vertex velocities
					pathline pl(*b, dtSim, dtParticle);
					pl.cell_to_vertex_interpolation(b, local_gcIds);

				});



				// do initial advection




				// balance workload

				master_kdt.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					b->particles = std::move(particles_hold);
				});
				




			
				// load based repartition
				bool wrap = false;
				size_t samples = 512;
				diy::kdtree(master_kdt, assigner, DIM, domain, &block::particles, samples, wrap);

				master_kdt.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					particles_hold = std::move(b->particles);
					dprint("particleshold %ld", particles_hold.size());
				});

				

				std::fill(gcIdToGid_local.begin(), gcIdToGid_local.end(), 0);
				std::fill(gcIdToGid_global.begin(), gcIdToGid_global.end(), 0);
				local_gcIds.clear();
				
				mpas1.velocityX.clear(); mpas1.velocityX.resize(mpas1.nCells*mpas1.nVertLevels);
				mpas1.velocityY.clear(); mpas1.velocityY.resize(mpas1.nCells*mpas1.nVertLevels);
				mpas1.velocityZ.clear(); mpas1.velocityZ.resize(mpas1.nCells*mpas1.nVertLevels);
				mpas1.vertVelocityTop.clear(); mpas1.vertVelocityTop.resize(mpas1.nCells*mpas1.nVertLevels);
				mpas1.zTop.clear(); mpas1.zTop.resize(mpas1.nCells*mpas1.nVertLevels);
				mpas1.zMid.clear(); mpas1.zMid.resize(mpas1.nCells*mpas1.nVertLevels);
				// dprint("mpas1.velocityX %ld %ld", mpas1.velocityX.size(), mpas1.nVertLevels);
				for (size_t i=0; i<particles_hold.size(); i++){
						if (particles_hold[i].predonly==2){
							
							int locCellIdx = particles_hold[i].glCellIdx;
							gcIdToGid_local[locCellIdx] = mpas1.gid;
							local_gcIds.push_back(locCellIdx);	
							// update the cell velocities in mpas1->velocityX; size nCells
							
							mpas1.velocityX.insert(mpas1.velocityX.begin()+locCellIdx*mpas1.nVertLevels, particles_hold[i].cell_data[0].begin(), particles_hold[i].cell_data[0].end());
							mpas1.velocityY.insert(mpas1.velocityY.begin()+locCellIdx*mpas1.nVertLevels, particles_hold[i].cell_data[1].begin(), particles_hold[i].cell_data[1].end());
							mpas1.velocityZ.insert(mpas1.velocityZ.begin()+locCellIdx*mpas1.nVertLevels, particles_hold[i].cell_data[2].begin(), particles_hold[i].cell_data[2].end());

							mpas1.vertVelocityTop.insert(mpas1.vertVelocityTop.begin()+locCellIdx*mpas1.nVertLevels, particles_hold[i].cell_data[3].begin(), particles_hold[i].cell_data[3].end());
							mpas1.zTop.insert(mpas1.zTop.begin()+locCellIdx*mpas1.nVertLevels, particles_hold[i].cell_data[4].begin(), particles_hold[i].cell_data[4].end());
							mpas1.zMid.insert(mpas1.zMid.begin()+locCellIdx*mpas1.nVertLevels, particles_hold[i].cell_data[5].begin(), particles_hold[i].cell_data[5].end());

						}
				}
				diy::mpi::all_reduce (world, gcIdToGid_local, gcIdToGid_global,  std::plus<int>());

				dprint("Finished prediction .............. ");

				block* block_ptr;
				master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					block_ptr = std::move(b);
				});

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

					b->gcIdToGid = std::move(gcIdToGid_global);
					std::string fname_graph = "graph.info";
					std::set<int> new_neighbors;

					b->create_links_from_gcIdToGid(fname_graph, new_neighbors);

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


				


			} // if (data_id == 15)

		} // looping through all in_data items
	} // outer while loop decaf get

	// if (world.rank()==0)
	// 	dprint("rank, %d, times trace lb gatherd, %f, %f, %f, ", world.size(), time_trace, time_lb, time_gatherd);

	// dprint("writing segments");
	// write segments to file in parallel
	// master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
	// dprint("writing");
	// b->parallel_write_segments(world, 0);
	// });
}

int main(int argc, char* argv[])
{
	
	dprint("In con");

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
				>> PosOption(seed_rate) >> PosOption(pred_percent)>> PosOption(skip_rate)))
    {
        if (world.rank() == 0)
        {
            fprintf(stderr, "Check ops usage \n", argv[0]);
            cout << ops;
        }
        return 1;
    }


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

	con(decaf, master, assigner, mpas1, pl, gp_file, dtSim, dtParticle);

	delete decaf;

	// MPI_Finalize(); // called by diy for consumers
}
