#include "streamlines.h"
#include "mpaso.h"
#include "flow.h"
#include "interpolators.h"
#include "nabo/nabo.h"
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

// generate seeds from the vertex list
void streamlines::generate_seeds(mpaso &mpasobj, int skipval, std::vector<double> &seeds_xyz){

    // std::cout<<mpasobj.xVertex[0]<<" "<<mpasobj.yVertex[0]<<" "<<mpasobj.zVertex[0]<<"\n";
    // std::cout<<mpasobj.zTopVertex[0]<<" "<<mpasobj.zTopVertex[0]<<"\n";


    std::vector<int> seed_level_ids = {10, 30, 50, 70, 90};

    int n_seed_verts = mpasobj.nVertices/skipval;
    int n_seed_levels = seed_level_ids.size();

    // std::cout<<mpasobj.nVertices/skipval;
    seeds_xyz.resize(3*n_seed_levels*n_seed_verts);


    double dc_b, dc_c, ratio; //dist from center of earth of bottom and current layer


    for (int i=0;i<n_seed_verts; i++){

        for (int j=0;j<n_seed_levels;j++){

            int vert_id = i*skipval;
            dc_c = radius + mpasobj.zTopVertex[mpasobj.nVertLevels*vert_id+seed_level_ids[j]];
            dc_b = radius + mpasobj.zTopVertex[mpasobj.nVertLevels*vert_id+0];
            ratio = dc_c/dc_b;

            seeds_xyz[i*n_seed_levels*3+j*3+0] = mpasobj.xVertex[vert_id] * ratio;
            seeds_xyz[i*n_seed_levels*3+j*3+1] = mpasobj.yVertex[vert_id] * ratio;
            seeds_xyz[i*n_seed_levels*3+j*3+2] = mpasobj.zVertex[vert_id] * ratio;

            //         std::cout<<vert_id<<" "<<j<<" "<<seeds_xyz.size()<<"\n";
        }

    }

}

streamlines::streamlines( mpaso &mpas1, std::vector<double> decaf_bar, std::string &op_file_name,
                         float init_t, float fin_t, float h, std::vector<double> &seeds_xyz){

    // read data







    double cx,cy,cz;
    const int K = 6;
    VectorXi nearest_idx(K);
    VectorXd dists2(K);
    VectorXd q(3);
    int num_iters = fin_t/h;

    mpas1.load_mesh_from_decaf_data_bar(decaf_bar);
    this->mpas1 = &mpas1;    
printf("nVertLevels in streamlines %ld\n", mpas1.nVertLevels);
    // for path_line path line this func needs to read in time as well. also need
    // time interval between writes to determine when fir_t and sec_t need to be updated
    // interpolation between times as well

    //mpas1.loadMeshFromNetCDF_CANGA(ip_file_name);

    // generate seed points put in list (cell centers ?)



    //    seeds_xyz.resize(n_seeds*3);
    //    seeds_xyz[0] = -5+5.61031e+06;
    //    seeds_xyz[1] = -5+237421;
    //    seeds_xyz[2] = -5+3.01006e+06;
    std::cout<<"starting\n";
    generate_seeds(mpas1, 10000, seeds_xyz);

    //    for (int i=0;i<seeds_xyz.size();i++){
    //        std::cout<<i<<" "<<seeds_xyz[i]<< "\n";
    //    }


    int n_seeds = seeds_xyz.size()/3;
    flow flow_lines(n_seeds);

    std::cout<<"nseeds "<<n_seeds<<"\n";
    // create kd trees

    MatrixXd C(3, mpas1.nCells);
    for (int i=0;i<mpas1.nCells;i++){
        C(0,i) = mpas1.xyzCell[i*3];
        C(1,i) = mpas1.xyzCell[i*3+1];
        C(2,i) = mpas1.xyzCell[i*3+2];
    }

    nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C);

    MatrixXd M(3, mpas1.nVertices);
    for (int i=0;i<mpas1.nVertices; i++){

        M(0,i) = mpas1.xVertex[i];
        M(1,i) = mpas1.yVertex[i];
        M(2,i) = mpas1.zVertex[i];

    }


    nns = Nabo::NNSearchD::createKDTreeLinearHeap(M);


    // stream generation line loop (handle out of boundary condition: if out, remove)

    bool terminate_stat;
    for (int j=0;j<n_seeds;j++){

        cx = seeds_xyz[j*3];
        cy = seeds_xyz[j*3+1];
        cz = seeds_xyz[j*3+2];

        flow_lines.fl_lengths[j] = num_iters;
        terminate_stat = 0;

        cur_t = init_t;
        // set line length = full
        flow_lines.fl_lengths[j] = num_iters;

	
        for (int i=0;i<num_iters;i++){

            q << cx, cy, cz;

            nns->knn(q, nearest_idx, dists2, K);

	

            // runge kutta to get next pos
            runge_kutta(cx, cy, cz, h, nearest_idx, terminate_stat);

            if (terminate_stat==1){
                flow_lines.fl_lengths[j] = i;
                break;
            }

            // if rk returns 1 set line length shorter and break

            flow_lines.flow_lines.push_back(cx);
            flow_lines.flow_lines.push_back(cy);
            flow_lines.flow_lines.push_back(cz);

            // increment time step and append pos to pathline
            cur_t += h;
//            std::cout<<i<<cx<<" "<<cy<<" "<<cz<<"\n";


        }



    }



    //    flow_lines.flow_lines.push_back(0); flow_lines.flow_lines.push_back(0); flow_lines.flow_lines.push_back(2);
    //    flow_lines.flow_lines.push_back(0); flow_lines.flow_lines.push_back(0); flow_lines.flow_lines.push_back(0);
    //    flow_lines.flow_lines.push_back(1); flow_lines.flow_lines.push_back(1); flow_lines.flow_lines.push_back(1);
    //    flow_lines.flow_lines.push_back(-2); flow_lines.flow_lines.push_back(3); flow_lines.flow_lines.push_back(5);

//    for (int i=0;i<flow_lines.fl_lengths.size();i++){
//        std::cout<<flow_lines.fl_lengths[i]<<" ";
//    }

    // write data to vtk
    flow_lines.write_flow_to_vtk(op_file_name);


}

void streamlines::runge_kutta(double &cx, double &cy, double &cz, float h, VectorXi &nearest_idx, bool &terminate_stat){

    //    std::vector<double> k1 = {0,0,0}, k2 = {0,0,0}, k3 = {0,0,0}, k4 = {0,0,0};
    //    std::vector<double> v1 = {0,0,0}, v2 = {0,0,0}, v3 = {0,0,0};

    Vector3d p = {cx, cy, cz};
    Vector3d pi = {cx, cy, cz};
    

//	std::cout<<"p0 p1 p2 c_vel nearest_idx"<<p[0]<<" "<< p[1]<<" "<<p[2]<<" "<<nearest_idx<<"\n";
    get_curpos_velocity_sl(p[0], p[1], p[2], c_vel, nearest_idx, terminate_stat);



    // get k_1
    Vector3d k1 = h*c_vel;

    // get v_1
    p = pi + (0.5 * k1);
    get_curpos_velocity_sl(p[0], p[1], p[2], c_vel, nearest_idx, terminate_stat);

    // get k_2
    Vector3d k2 = h*c_vel;

    // get v_2
    p = pi + (0.5 * k2);
    get_curpos_velocity_sl(p[0], p[1], p[2], c_vel, nearest_idx, terminate_stat);

    // get k_3
    Vector3d k3 = h*c_vel;

    //get v_3
    p = pi + k3;
    get_curpos_velocity_sl(p[0], p[1], p[2], c_vel, nearest_idx, terminate_stat);


    // get k_4
    Vector3d k4 = h*c_vel;

    // get p_next
    Vector3d p_next = pi + (k1 + k2 + k3 + k4)/6;

    cx = p_next[0];
    cy = p_next[1];
    cz = p_next[2];



}
void streamlines::get_curpos_velocity_sl(double cx, double cy, double cz, Vector3d &c_vel,
                                         VectorXi &nearest_idx, bool &terminate_stat){

    // get depth from cx, cy, cz
    double depth = std::sqrt(cx*cx + cy*cy + cz*cz) - radius;

    //    std::cout<<radius<<"\n";

    //    double dx = xyzCell[600*3],dy = xyzCell[600*3+1], dz = xyzCell[600*3+2];
    //    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    //    std::cout<<dist<<"\n";
    //    std::cout<<dist - radius;

    VectorXi nearest_cell_idx(1);
    VectorXd dists2(1);
    VectorXd q(3);

    q << cx, cy, cz;
    // interpolate vertically and get positions and velocity (in values)
    std::vector<double> values;
    values.resize(6*nearest_idx.size());

    // create another kd tree for cells. use to check for alt

    nns_cells->knn(q, nearest_cell_idx,dists2, 1);
//std::cout<<"reaching here 1 "<<mpas1->nVertLevels;//*nearest_cell_idx[0]+0;


    double top_boun = mpas1->zTop[mpas1->nVertLevels*nearest_cell_idx[0]+0];

    double botom_boun = mpas1->zTop[mpas1->nVertLevels*nearest_cell_idx[0]+99];



    if (depth> top_boun || depth< botom_boun){
        // terminate flow line
        terminate_stat = 1;
//        std::cout<<top_boun<<" "<<depth<<" "<<botom_boun<<" "<<"ya"<<"\n";
//        std::cout<<"terminating line\n";

    }else{
//        std::cout<<"all good\n";
        interpolate_vertically(mpas1->nVertLevels, mpas1->zTopVertex, nearest_idx, values, depth,
                               mpas1->xVertex, mpas1->yVertex, mpas1->zVertex,
                               mpas1->velocityXv, mpas1->velocityYv, mpas1->velocityZv);


    // get nearest vertex, if any neiboring cell zero, then horizontal boundary

    // interpolate horizontally
    interpolate_horizontally(cx, cy, cz, values, c_vel);

   }


}
void streamlines::compute(){




}

void streamlines::write(){


}

streamlines::~streamlines(){

}
