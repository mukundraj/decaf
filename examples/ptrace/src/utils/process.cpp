#include "process.h"
#include "block.h"
#include "streamline.h"
#include "misc.h"



void remote_enq_seeds(
        block* b,
        const diy::Master::ProxyWithLink&   cp,
        const diy::Assigner&                assigner
){


    std::map<diy::BlockID, std::vector<EndPt>>   outgoing_seeds;
    for (int i=0; i<b->particles.size(); i++){

        EndPt &endpt = b->particles[i];

        int dest_gid = b->gcIdToGid_init[endpt.glCellIdx+1];


        int dest_proc = assigner.rank(dest_gid);
        diy::BlockID dest_block = {dest_gid, dest_proc};

        outgoing_seeds[dest_block].push_back(endpt);
    }

    // enqueue into cp
    for (auto seed: outgoing_seeds){
		cp.enqueue(seed.first, seed.second);
	}

	
    
}

void remote_deq_seeds(block* b, const diy::Master::ProxyWithLink& cp){



    b->particles.clear();

    vector<int> in;
	cp.incoming(in);
	for (int i = 0; i < in.size(); i++)
	{
		// if (cp.incoming(in[i]).buffer.size() > 0)
		if (cp.incoming(in[i]).size() > 0)
		{
			std::vector<EndPt> incoming_endpts;
			cp.dequeue(in[i], incoming_endpts);
			// b->process_halo_req(incoming_halo, framenum);
			// dprint("cp.gid %d, numdeq %d", cp.gid(), incoming_endpts.size());

			for (size_t j=0; j<incoming_endpts.size(); j++){
                b->particles.push_back(incoming_endpts[j]);
				// b->copy_endpt_to_cell_data(incoming_endpts[j]);
			}
		}
	}

}

void send_sampled_seeds_to_origin(diy::Master &master, const diy::Assigner&                assigner){

    master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
					
        remote_enq_seeds(b, cp, assigner);

		

	});	


    bool remote = true;
    master.exchange(remote);
    master.foreach(&remote_deq_seeds);

}


void remote_enq_data(
        block* b,
        const diy::Master::ProxyWithLink&   cp,
        const diy::Assigner&                assigner){

    std::map<diy::BlockID, std::vector<EndPt>>   outgoing_data;

	// dprint("enqueing to sort b->local_gcIds %ld, gid %d", b->local_gcIds_init.size(), cp.gid());

    for (size_t i=0; i<b->local_gcIds_init.size(); i++){
       
	   int gcId = b->local_gcIds_init[i];
	   if (b->gcIdToGid[gcId] != cp.gid()){

			int cell_idx = gcId - 1;
			EndPt cell_cen;
			cell_cen[0] = b->xCell[cell_idx];
			cell_cen[1] = b->yCell[cell_idx];
			cell_cen[2] = b->zCell[cell_idx];
			cell_cen.predonly = 2; // indicates this is the cell center point
			cell_cen.glCellIdx = cell_idx;
			cell_cen.source_gid = cp.gid();

			b->copy_cell_data_to_endpt(cell_cen, cell_cen.glCellIdx);

			int dest_gid = b->gcIdToGid[cell_cen.glCellIdx+1];
			int dest_proc = assigner.rank(dest_gid);
			diy::BlockID dest_block = {dest_gid, dest_proc};
			outgoing_data[dest_block].push_back(cell_cen);
			// dprint("cpgid %d dest_gid %d", cp.gid(), dest_gid);

	   }
        
    }

    // enqueue into cp
    for (auto celldata: outgoing_data){
		// dprint("from %d to %d : %ld", cp.gid(), celldata.first.dest_gid, outgoing_data.size());
		cp.enqueue(celldata.first, celldata.second);
	}

}

void remote_deq_data(block* b, const diy::Master::ProxyWithLink& cp){

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
                // b->particles.push_back(incoming_endpts[j]);
				if (incoming_endpts[j].glCellIdx == 5550 - 1)
					dprint("found in %d", cp.gid());
				b->copy_endpt_to_cell_data(incoming_endpts[j]);
			}
		}
	}

}




void ghost_exchange(diy::Master &master, const diy::Assigner& assigner, bool pred_frame){

    master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {	

        // b->gcIdToGid = std::move(gcIdToGid_global);
        enqueue_ghost_cells_to_nbrs(b, cp, assigner, pred_frame);
        
    });

    master.exchange();

    // if (world.rank()==0)
    //     dprint("after exchange, framenum %d, local_gcIds %ld", framenum, local_gcIds.size());	


    master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {


        deque_ghost_cells_from_nbrs(b, cp);

        // // compute vertex velocities
        // streamline sl(*b, dtSim, dtParticle);
        // sl.cell_to_vertex_interpolation(b, local_gcIds);

    });

}

void enqueue_ghost_cells_to_nbrs(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int pred_frame){

	std::vector<int> &local_gcids = (pred_frame==true) ? b->local_gcIds_init : b->local_gcIds; 
	std::vector<int> &gcidtogid = (pred_frame==true) ? b->gcIdToGid_init : b->gcIdToGid; 

	diy::Link*    l = cp.link();


	
	std::map<diy::BlockID, std::vector<EndPt>>   outgoing_cells;

	// dprint("before enqueing ghost b->local_gcIds_init %ld, b->indexToCellID %ld, l", b->local_gcIds_init.size(), b->indexToCellID.size());
	for (size_t i = 1; i < local_gcids.size(); i++)
	{	
		int gcId = local_gcids[i];
		// if (i<10)
			// dprint("i %d size %ld, gid %d", i, b->cell_nbrs[gcId].size(), cp.gid());
			// dprint("rank %d gcId %d", cp.gid(), gcId);

		
		for (size_t j=0; j<b->cell_nbrs[gcId].size();j++){

			int nbr_gcId = b->cell_nbrs[gcId][j]; 

			// dprint("rank %d nbr_gcId %d, b->gcIdtoGid %ld", cp.gid(), nbr_gcId, b->gcIdToGid.size());
			int nbr_gid = gcidtogid[nbr_gcId];

			int link_id = l->find(nbr_gid);

			if (link_id > 0){

				diy::BlockID dest_block = l->target(link_id);

				EndPt endpt;
				endpt.glCellIdx = gcId - 1;

				// dprint("b->velocityX %ld", b->velocityX.size());
				b->copy_cell_data_to_endpt(endpt, endpt.glCellIdx);
				outgoing_cells[dest_block].push_back(endpt);
			}

			
			// if (nbr_gid != b->gid){
			// 	// send info of gcId to nbr_gid, store in cell_data
			// 	EndPt endpt;
			// 	endpt.glCellIdx = gcId - 1;

				
			// 	// dprint("b->velocityX %ld", b->velocityX.size());
			// 	b->copy_cell_data_to_endpt(endpt, endpt.glCellIdx);
				
			// 	// if (b->framenum >2)
			// 	// continue;

			// 	int dest_gid = nbr_gid;
			// 	// if (dest_gid<0 || dest_gid>32)
			// 	// 	dprint("alert rank %d, dest_gid %d", cp.gid(), dest_gid);
			// 	int dest_proc = assigner.rank(dest_gid);
			// 	diy::BlockID dest_block = {dest_gid, dest_proc};

			// 	outgoing_cells[dest_block].push_back(endpt);

				
				
			// }
		}
		// if (i > local_gcIds.size()-2)
		// 	dprint("rank %d gcId %d i %ld local_gcIds %ld", cp.gid(), gcId, i, local_gcIds.size());

	}
	// dprint("rank %d done", cp.gid());
	// else if not in rebalance interval
	// loop over keys of map<int, std::vector<int>> epoch_nbrs and cread a endpt for each nbr

	// dprint("before enqueing");

	for (auto elem: outgoing_cells){
		cp.enqueue(elem.first, elem.second);
	}

	// dprint("after enqueing gid %d", cp.gid());
	
	// dprint("gid %d outsize %ld", b->gid, outgoing_cells.size());
	// if (b->gid==2){
	// 	for (auto i : outgoing_cells){
	// 		dprint("outs %d", i.second.size());
	// 	}

	// }
}

void deque_ghost_cells_from_nbrs(block *b, const diy::Master::ProxyWithLink &cp){

	// vector<int> in;
	// cp.incoming(in);
	// for (int i = 0; i < in.size(); i++)
	// {
	// 	if (cp.incoming(in[i]).buffer.size() > 0)
	// 	{
	// 		std::vector<EndPt> incoming_endpts;
	// 		cp.dequeue(in[i], incoming_endpts);

	// 		for (size_t j=0; j<incoming_endpts.size(); j++){
	// 			b->copy_endpt_to_cell_data(incoming_endpts[j]);
	// 		}
	// 	}
	// }

	diy::Link* l = cp.link();
	for (int i = 0; i < l->size(); ++i)
    {
        int gid = l->target(i).gid;
        if (cp.incoming(gid).size())
        {
            std::vector<EndPt> incoming_endpts;
			cp.dequeue(cp.link()->target(i).gid, incoming_endpts);

			for (size_t j=0; j<incoming_endpts.size(); j++){
				b->copy_endpt_to_cell_data(incoming_endpts[j]);
			}
        }
    }
	
	
}


void interpolate_cell_to_vertex(diy::Master &master, bool pred_frame){

    master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
        streamline sl;
		
	    sl.cell_to_vertex_interpolation(b, pred_frame);
	});	
    
}

// sorting cell data after rebalance to current partitions
void sort_data(diy::Master &master, const diy::Assigner& assigner){

	

    master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {

		b->move_data_to_regular_position(cp.gid());

        remote_enq_data(b, cp, assigner);

		
	});	

    bool remote = true;
    master.exchange(remote);
    master.foreach(&remote_deq_data);

}

void print_global_num_particles(std::string &str, diy::mpi::communicator& world, diy::Master &master, int framenum){

	size_t total_remaining_particles = 0;



	master.foreach ([&](block *b, const diy::Master::ProxyWithLink &cp) {
		diy::mpi::reduce(world, b->particles.size(), total_remaining_particles, 0, std::plus<size_t>());
	});	

	world.barrier();
	if (world.rank()==0)
		dprint("%s %ld, frame %d", str.c_str(), total_remaining_particles, framenum);

}