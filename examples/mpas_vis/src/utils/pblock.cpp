
#include "pblock.h"
#include <fstream>
#include <sstream>
#include "misc.h"
#include "nabo/nabo.h"
#include <Eigen/Dense>
#include <diy/mpi.hpp>

EndPt::EndPt()
        {
            pid      = 0;
            sid      = 0;
	    step     = 0; // initially all particles are at step zero
	    gpid     = 0; // need to set this based on the block of initialization
        }
    

    // following constructor defined out of line because references Segment, which needed
// to be defined first
EndPt::
EndPt(Segment& s)                       // extract the end point of a segment.  CHECK also pass start_step here?
{
    pid = s.pid;
    sid = s.sid;
    step = s.start_step + s.pts.size();
    gpid = s.gpid;
    pt.coords[0] = s.pts.back().coords[0];
    pt.coords[1] = s.pts.back().coords[1];
    pt.coords[2] = s.pts.back().coords[2];
    pt.coords[3] = s.pts.back().coords[3];
}



    Segment::Segment()
        {
            pid      = 0;
            sid      = 0;
        }
    Segment::Segment(EndPt& p  )       // construct a segment from one point.
        {
            pid      = p.pid;
            sid      = p.sid;
            Pt pt    = { p[0], p[1], p[2] };
	    start_step = p.step;
            pts.push_back(pt);
	    gpid = p.gpid;
        }

    // // whether end point is inside given bounds
    // bool Segment::inside(const int lb[3], const int ub[3]) const
    //     {
    //         for (int i = 0; i < 3; i++)
    //             if (pts.back().coords[i] < lb[i] || pts.back().coords[i] >= ub[i] - 1)
    //                 return false;
    //         return true;
    //     }

void PBlock::read_cell_g_neighbors(){

    std::string ip_file = "graph.info";

    std::ifstream file(ip_file.c_str());
    // file.clear();
    std::string line;
    int temp;
    std::getline(file,line); // reading the header line
    std::istringstream iss(line);
    iss >>  temp;
    cell_g_neighbors.resize(temp+1); 


    int cur_cell_1idx = 1;
    while (std::getline(file,line)) { 
        std::istringstream iss2(line);
        while(iss2 >>  temp){

            cell_g_neighbors[cur_cell_1idx].push_back(temp);
        }
        cur_cell_1idx++;
        
        // fprintf(stderr, " cur_cell_1idx %d\n", cur_cell_1idx);
        
    }

}

void PBlock::generate_domain_decomposition_graph(std::string &filename, int nblocks){

    std::string ip_file = "graph.info.part."+itos(nblocks);
    int ctr=0, temp;
    std::ifstream file;
    file.open(ip_file);
    std::string line;
    ctr=1;
    while (std::getline(file,line)) { 
        std::istringstream iss(line);
        iss >>  temp;
        cgid_to_bid[ctr] = temp;

        ctr++;
        
    }
    file.close();
    

}

void PBlock::compute_cellIndex(int cblock, int nblocks){

    // now populate the cellID_to_bid vector
    std::string ip_file = "graph.info.part."+itos(nblocks);
    std::ifstream file;
    file.clear();
    
    int g_idx=0, cblock_idx=0, temp;
    file.open(ip_file);

    std::string line;

    while (std::getline(file,line)) { 
        std::istringstream iss(line);
        iss >>  temp;
        // fprintf(stderr, " temp %d %d\n", temp, cblock);
        if (temp==cblock){
            cellIndex[g_idx+1] = cblock_idx;
            cblock_idx++;
        }


        g_idx++; 
    }

    file.close();


}

int PBlock::get_bid_for_pt(double *coords){
    
    int bid;    
    
    // find nearest cell
    Eigen::VectorXi nearest_cell_idx(1);
    Eigen::VectorXd dists2(1);
    Eigen::VectorXd q(3);

    q << coords[0], coords[1] , coords[2];      
    nns_cells->knn(q, nearest_cell_idx,dists2, 1);
    // dprint("nearest_cell_idx cgid %d %d %d %d", nearest_cell_idx[0], indexToCellID[nearest_cell_idx[0]], cgid_to_bid[indexToCellID[nearest_cell_idx[0]]], cgid_to_bid[4463]);
    // find the bid of the cell
    bid = cgid_to_bid[indexToCellID[nearest_cell_idx[0]]];
    return bid;     
}

void PBlock::initialize_seeds(int skipval, diy::mpi::communicator &world, vector<EndPt> &carryover_particles){

        double dc_t, dc_c, ratio, dc_b, dc_base; //dist from center of earth of top and current layer


            // int skipval = 5000;
        // int skipval = 2500;
        // std::vector<int> seed_level_ids = {10, 30, 50, 70, 90};
        std::vector<int> seed_level_ids = { 10};
        int n_seed_verts = nVertices/skipval;
        int n_seed_levels = seed_level_ids.size();

        int n_particles = n_seed_verts * n_seed_levels;
        diy::mpi::all_gather(world, n_particles, global_nP);
        global_start_ids.resize(global_nP.size());
        // fprintf(stderr, "global np size %ld \n", global_nP.size());
        
        int n_particles_global = global_nP[0];
        for (int k = 1; k<global_nP.size(); k++){
            global_start_ids[k] = global_start_ids[k-1]+global_nP[k-1];
            // fprintf(stderr, "global_start_pos %d %d\n", k, global_start_pos[k]);
            n_particles_global += global_nP[k];
        }
        global_trace_sizes.resize(n_particles_global);
        dprint("n_seed_verts %d gid %d", n_seed_verts, gid);
        for (int i=0;i<n_seed_verts; i++){

            for (int j=0;j<n_seed_levels;j++){

                int vert_id = i*skipval;
                // if (gid==1){
                //  dprint("gid 1 seed_level_ids[j] %d zTopVertex %f",seed_level_ids[j],mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+seed_level_ids[j]]);
                // }


                dc_c = radius + zTopVertex[nVertLevels*vert_id+seed_level_ids[j]];
                dc_t = radius + zTopVertex[nVertLevels*vert_id+0];
                dc_b = radius + zTopVertex[nVertLevels*vert_id+99];




                dc_base = std::sqrt(xVertex[vert_id]*xVertex[vert_id] + yVertex[vert_id]*yVertex[vert_id] + zVertex[vert_id]*zVertex[vert_id]);
                ratio = dc_c/dc_base;

                EndPt p;
                p.pid = init;
                p.sid = init;
                p.step = 0;
                p.gpid = global_start_ids[gid] + p.pid;
                p[0] = xVertex[vert_id] * ratio;
                p[1] = yVertex[vert_id] * ratio;
                p[2] = zVertex[vert_id] * ratio;
                p[3] = 0; // later p.step * h


                if (gid==1){

                    dc_base = std::sqrt(-5825858.330663*(-5825858.330663) + -282189.163605*(-282189.163605) + -2563650.318867*(-2563650.318867));
                    ratio = dc_c/dc_base;
                    p[0] = -5825858.330663 * ratio;
                    p[1] = -282189.163605 * ratio;
                    p[2] = -2563650.318867 * ratio;
                    p[3] = 0; // later p.step * h
                }


                // dprint("vlid gvid bid %d %d %d", vert_id, indexToVertexID[vert_id], gid);

                carryover_particles.push_back(p);
                // init++; // moved to consumer loop



            }
        }

}


void* create_block()                      { return new PBlock;  }
void  destroy_block(void* b)              { delete static_cast<PBlock*>(b);  }
void  save_block(const void* b,
                 diy::BinaryBuffer& bb)   { 
						diy::save(bb, *static_cast<const PBlock*>(b));
											
						  
					}
void  load_block(void* b,
                 diy::BinaryBuffer& bb)   { diy::load(bb, *static_cast<PBlock*>(b));  }