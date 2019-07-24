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


void enq_halo_info(block *b, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner){


		Halo h;

		std::vector<Halo> halos = b->get_halo_info();
		for (size_t i=0; i<halos.size(); i++){

			int dest_gid            = halos[i].dest_gid;
			int dest_proc           = assigner.rank(dest_gid);
			diy::BlockID dest_block = {dest_gid, dest_proc};
			cp.enqueue(dest_block, halos[i]);

		}
		
}


void process_halo_req(block *b, const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, int framenum){


		vector<int> in;
		cp.incoming(in);
		for (int i = 0; i < in.size(); i++)
		{	
			if (cp.incoming(in[i]).buffer.size() > 0)
			{
				Halo incoming_halo;
				cp.dequeue(in[i], incoming_halo);
				b->process_halo_req(incoming_halo, framenum);

				int dest_gid            = incoming_halo.src_gid;
				int dest_proc           = assigner.rank(dest_gid);
				diy::BlockID dest_block = {dest_gid, dest_proc};
				cp.enqueue(dest_block, incoming_halo);

			}
		}


}

void deq_halo_info(block *b, const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, int framenum){

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
void  process_halo_dynamic(block *b, const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, int framenum){

	// b->process_halo_dynamic(incoming_halo, framenum);

}

void deq_halo_dynamic(block *b, const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, int framenum){

	// b->update_halo_dynamic(incoming_halo, framenum);

}

// consumer
void con(Decaf* decaf, diy::Master &master, diy::RoundRobinAssigner &assigner, block& mpas1, pathline &pl)
{   


	

	// mpaso mpas_c;
	vector< pConstructData > in_data;
	diy::mpi::communicator world(decaf->con_comm_handle());

	// recomputed particle file ("particles2.nc") includes global cell indices for seeds
	std::string fname_particles = "particles2.nc", fname_graphinfo = "graph.info.part."+std::to_string(world.size());	

	while (decaf->get(in_data))
	{
		int n_velocityX = 0;
		std::vector<int> indexToCellID;
		std::vector<double> data_bar_dbl;
		std::vector<int> data_bar_int;
		int fir_cell_idx = 999;
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
			mpas1.update_data(data_id, framenum, data_bar_int, data_bar_dbl);

			if (framenum==1 && data_id==15){

				// create links
				diy::Link* link = new diy::Link; 

				std::set<int> links;
				mpas1.create_links_mpas(fname_graphinfo, links, world);

				for (auto link_gid : links) {
					diy::BlockID  neighbor;
					neighbor.gid  = link_gid;                     // gid of the neighbor block
					neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
					link->add_neighbor(neighbor);                // add the neighbor block to the link
				}

				// add links and block to master
				master.add(mpas1.gid, &mpas1, link);

				master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
                {	
					enq_halo_info(b, cp, assigner);
				});
				master.exchange();
				master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
                {	
					process_halo_req(b, cp, assigner, framenum);
				});
				master.exchange();
				master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
                {	
					deq_halo_info(b, cp, assigner, framenum);
				});

				// initialize seeds
				// mpas1.generate_new_particle_file();
				mpas1.init_seeds_mpas(fname_particles);
			
			}else if (framenum>1 && data_id==13){
				
				// dynamic halo exchange
				master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
                {	
					process_halo_dynamic(b, cp, assigner, framenum);
				});
				master.exchange();
				master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
                {	
					deq_halo_dynamic(b, cp, assigner, framenum);
				});

				//do{
				master.foreach([&](block* b, const diy::Master::ProxyWithLink& cp)
                {		
                		// update field
                		pl.update_velocity_vectors(*b, framenum);

                		// trace particles
                    	pl.compute_epoch(b);
				});
				// master.exchange();
				//while(reduce of all remaining particles not zero)

			}


			
			

			
		}
	}

}

int main(){

	

	block mpas1;	
	double dtSim=7200, dtParticle=300;
	pathline pl(mpas1, dtSim, dtParticle);

		// define the workflow
	Workflow workflow;
	//make_wflow(workflow);
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");



	MPI_Init(NULL, NULL);


	// create decaf
	Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

	diy::mpi::communicator world(decaf->con_comm_handle());
		// dprint("diy world size: %d", world.size() );

	// if (world.rank() == 0)
 //    {
        
 //        volatile  int dwait=0;
 //        fprintf(stderr , "pid %ld  waiting  for  debugger\n"
 //            , (long)getpid ());
 //            while(dwait==0) { /*  change  ’i’ in the  debugger  */ }
 //    }
 //    world.barrier();

	diy::FileStorage storage("./DIY.XXXXXX");
	int                       nblocks    = world.size();
	int                       threads    = 1;
	int                       mem_blocks = -1;

	diy::Master               master(world,
		threads,
		mem_blocks,
		NULL,
		NULL,
		&storage,
		NULL,
		NULL);

	diy::RoundRobinAssigner assigner(world.size(), nblocks);

	std::vector<int> gids;                     // global ids of local blocks
	assigner.local_gids(world.rank(), gids); // get the gids of local blocks
	mpas1.gid = gids[0];


	con(decaf, master, assigner, mpas1, pl);


	delete decaf;

	// MPI_Finalize(); // called by diy for consumers

}


// switch(d_metadata.getData()){

			// 	case 0: fprintf(stderr, "Recv xCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 1: fprintf(stderr, "Recv yCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 2:	fprintf(stderr, "Recv zCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 3:	fprintf(stderr, "Recv xVertex %d,\n", d_metadata.getData());
			// 			break;
						
			// 	case 4:	fprintf(stderr, "Recv yVertex %d,\n", d_metadata.getData());
			// 			break;

			// 	case 5:	fprintf(stderr, "Recv zVertex %d,\n", d_metadata.getData());
			// 			break;

			// 	case 6:	fprintf(stderr, "Recv indexToVertexID %d,\n", d_metadata.getData());
			// 			break;

			// 	case 7:	fprintf(stderr, "Recv indexToCellID %d,\n", d_metadata.getData());
			// 			break;

			// 	case 8:	fprintf(stderr, "Recv verticesOnEdge %d,\n", d_metadata.getData());
			// 			break;

			// 	case 9:	fprintf(stderr, "Recv cellsOnVertex %d,\n", d_metadata.getData());
			// 			break;

			// 	case 10:	fprintf(stderr, "Recv verticesOnCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 11:	fprintf(stderr, "Recv velocityXv %d, fno %d\n", d_metadata.getData(), frame_num.getData());
			// 			break;

			// 	case 12:	fprintf(stderr, "Recv velocityYv %d, fno %d\n", d_metadata.getData(), frame_num.getData());
			// 			break;

			// 	case 13:	fprintf(stderr, "Recv velocityZv %d,\n", d_metadata.getData());
			// 			break;

			// 	case 14:	fprintf(stderr, "Recv nEdgesOnCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 15:	fprintf(stderr, "Recv maxLevelCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 16:	fprintf(stderr, "Recv boundaryVertex %d,\n", d_metadata.getData());
			// 			break;

			// 	case 17:	fprintf(stderr, "Recv vertVelocityTop %d,\n", d_metadata.getData());
			// 			break;

			// 	case 18:	fprintf(stderr, "Recv cellsOnCell %d,\n", d_metadata.getData());
			// 			break;

			// 	case 19:	fprintf(stderr, "Recv zTop %d,\n", d_metadata.getData());
			// 			break;

			// 	case 20:	fprintf(stderr, "Recv zMid %d,\n", d_metadata.getData());
			// 			break;

			// 	default:fprintf(stderr, "Exiting! %d,\n", d_metadata.getData()); 
			// 			exit(0);
			// 			break;


			// }

			// if(d_metadata.getData()==19)
			// {	
			// 	VectorFliedd d_data_bar = in_data[i]->getFieldData<VectorFliedd>("data_d");

			// 	data_bar = d_data_bar.getVector();
			// 	printf("Incomingg container %f metadata2 %d\n", data_bar[data_bar.size()-1], d_metadata.getData());

				
			// }else if (d_metadata.getData()==7){
			// 	VectorFieldi d_data_bar = in_data[i]->getFieldData<VectorFieldi>("data_i");

			// 	data_bar_int = d_data_bar.getVector();
			// 	printf("Incomingg container %d metadata2 %d\n", data_bar_int[data_bar.size()-1], d_metadata.getData());

			// }else if (d_metadata.getData()==0){

			// 	VectorFliedd d_data_bar = in_data[i]->getFieldData<VectorFliedd>("data_d");
			// 	data_bar = d_data_bar.getVector();
			// 	printf("Incomingg container %f metadata2 %d\n", data_bar[data_bar.size()-1], d_metadata.getData());
			// }