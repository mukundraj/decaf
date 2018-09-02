#include "interpolators.h"

// *** need to linearly interpolate positions( x, y and z) and values ***

void interpolate_vertically(size_t nVertLevels, std::vector<double> &zTopVertex, Eigen::VectorXi &nearest_idx,
                            std::vector<double> &values, double depth,
                            std::vector<double> &xVertex, std::vector<double> &yVertex, std::vector<double> &zVertex,
                            std::vector<double> &velocityXv, std::vector<double> &velocityYv, std::vector<double> &velocityZv){

    const double radius = 6371220.;

    int vert_id, bot_idx, top_idx;

    double dc_b, dc_m, ratio; //distance from center for bottom/mid layer.

    //    for (int i=0;i<nearest_idx.size();i++){
    //        std::cout<<nearest_idx[i]<<"\n";
    //    }

    //    for (int i=0; i<nVertLevels; i++){
    //        std::cout<<zTopVertex[nVertLevels*3700+i]<<"\n";
    //    }

    int nnsize = nearest_idx.size();
    // for each id



    for (int i=0;i<nnsize;i++){
        vert_id = nearest_idx[i];

        //get upper and lower
        bot_idx = bin_search_index(zTopVertex, nVertLevels*vert_id,nVertLevels*vert_id+nVertLevels, depth);// pass depth
        bot_idx = bot_idx % nVertLevels;
        top_idx = bot_idx - 1;

        //        std::cout<<bot_idx<<" "<<top_idx<<"\n";

        // interpolate x, y, z, Vx, Vy, Vz

        //        store interpolated values in values
        values[6*i] = linear_inter(depth, zTopVertex[nVertLevels*vert_id+bot_idx], zTopVertex[nVertLevels*vert_id+top_idx],
                velocityXv[vert_id*nVertLevels+bot_idx],velocityXv[vert_id*nVertLevels+top_idx]);
        values[6*i+1] = linear_inter(depth, zTopVertex[nVertLevels*vert_id+bot_idx], zTopVertex[nVertLevels*vert_id+top_idx],
                velocityYv[vert_id*nVertLevels+bot_idx],velocityYv[vert_id*nVertLevels+top_idx]);
        values[6*i+2] = linear_inter(depth, zTopVertex[nVertLevels*vert_id+bot_idx], zTopVertex[nVertLevels*vert_id+top_idx],
                velocityZv[vert_id*nVertLevels+bot_idx],velocityZv[vert_id*nVertLevels+top_idx]);

        dc_b = radius + zTopVertex[nVertLevels*vert_id+bot_idx];
        dc_m = radius + depth;
        ratio = dc_m/dc_b;

        values[6*i+3] = xVertex[vert_id] * ratio;
        values[6*i+4] = yVertex[vert_id] * ratio;
        values[6*i+5] = zVertex[vert_id] * ratio;

        //        std::cout<<dc_m<<" "<<dc_b<<" "<<dc_m/dc_b<<"\n";


    }
    //    for (int i=0;i<values.size();i++){
    //        std::cout<<i<<" "<<values[i]<<"\n";
    //    }
}

double linear_inter(double x, double x1, double x2, double q00, double q01){
    return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
}

void interpolate_horizontally(double cx, double cy, double cz,
                              std::vector<double> &values, Eigen::Vector3d &c_vel){


    size_t n = values.size()/6;

    std::vector<double> b(6);
    std::vector<Eigen::Vector3d> s(n);
    std::vector<Eigen::Vector3d> e(n);
    std::vector<Eigen::Vector3d> v(n);
    std::vector<Eigen::Vector3d> vel(n);
    Eigen::Vector3d p(cx, cy, cz);


    for (int i=0;i<n;i++){
        vel[i][0] = values[6*i+0];
        vel[i][1] = values[6*i+1];
        vel[i][2] = values[6*i+2];

        v[i][0] = values[6*i+3];
        v[i][1] = values[6*i+4];
        v[i][2] = values[6*i+5];
    }

    for (int i = 0; i < n; i++) {

        const size_t ip = (i + 1) % n;

        s[i] = v[i] - p;
        e[i] = v[ip] - v[i];
    }

    std::vector<double> A(n, 1.0);
    std::vector<double> C(n);

    for (size_t i=0;i<n;i++){
        const size_t im = (i + n - 1) % n;
        for (size_t j=0;j<n;j++){
            if (j != im && j != i) {

                const size_t jp = (j + 1) % n;
                //                        std::cout<<s[j].cross(s[jp]).norm();

                A[i] *= 0.5 * s[j].cross(s[jp]).norm();
            }

        }
        C[i] = 0.5 * e[i].cross(-e[im]).norm();

    }
    std::vector<double> w(n);
    double W = 0.0;

    for (size_t i = 0; i < n; ++i) {

        w[i] = C[i] * A[i];
        W += w[i];
    }
    if (fabs(W)>=0){
        assert(fabs(W) > 0.0);
        const double invW = 1.0 / W;

        for (size_t i = 0; i < n; ++i){
            b[i] = w[i] * invW;
        }

        // remember to fill in c_vel
        c_vel<<0,0,0;
        for (size_t i = 0; i < n; ++i){
            c_vel += b[i]*vel[i];
        }
    }else{
        c_vel<<0,0,0;
    }

    //    std::cout<<c_vel;


}



int bin_search_index(std::vector<double> &zTopVertex, int lo, int hi, double val){
    int mid, next=hi;

    while (lo < hi){
        mid = lo + (hi-lo) / 2;

        if (zTopVertex[mid]==val){
            return mid+1;
        } else if (zTopVertex[mid]>val && zTopVertex[mid+1]<val){

            return mid+1;

        } else if (zTopVertex[mid-1]>val && zTopVertex[mid]<val){

            return mid;

        }else if (val < zTopVertex[mid]){
            lo = mid+1;
            next = lo;
        }else{
            next = mid;
            hi = mid-1;

        }

    }
    return next;
}


