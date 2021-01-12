
// Prepares a dense seed particles file for LIC generating LIC images by essentially identifying the cells
// that each source/high res seed is part of in the low/target resolution test case.
// Input: 
//      - A particles_lic.nc file with dense seeding generating using make_particles python script for a high resolution test case


#include "particles_io.h"
#include <string>
#include <vector>
#include "misc.h"



bool point_in_cell(int nVertices, double *xv, double *yv, double *zv, double xp, double yp, double zp){

    

    return true;

}

void get_cell_ids_for_seeds(seed_data &data){

    for (size_t i=0; i<data.xPos.size(); i++){

        int cell_id;

        //point_in_cell

        break;

    }


}


int main(){

    // params 
    std::string ip_particles = "/home/mraj/tests/baseline2/mpas-decaf-ptrace/MPAS-Model/testing_and_setup/MPAS-O_V6.0_EC60to30/particles/particles32.nc";

    std::string ip_mesh = "/home/mraj/tests/baseline2/mpas-decaf-ptrace/MPAS-Model/testing_and_setup/MPAS-O_V6.0_QU240_decaf_flowvis/output.nc";


    // local vars
    seed_data data;

    
    // read source particles.nc 
    load_particle_data(ip_particles, data);

    // read cell center locations of target mesh
    load_MPAS_cells_coords(ip_mesh, data);

    // loop over each particle and identify which target mesh cell the particle belongs to
    get_cell_ids_for_seeds(data);

    dprint("num_particles %ld", data.xPos.size());

    // write out particles_lic.nc file with updated cell ids along with other particle information


    return 0;
}



