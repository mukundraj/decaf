#ifndef PTRACE_H
#define PTRACE_H

//---------------------------------------------------------------------------
//
// diy2-vtk7 parallel particle advection trace classes
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------
// #include <diy/mpi.hpp>
// #include <diy/master.hpp>
// #include <diy/reduce-operations.hpp>
// #include <diy/decomposition.hpp>
// #include <diy/mpi/datatypes.hpp>
// #include <diy/io/bov.hpp>
// #include <diy/pick.hpp>
// #include <diy/reduce.hpp>
// #include <diy/partners/merge.hpp>
#include <diy/point.hpp>
#include <diy/serialization.hpp>

using namespace std;

// typedef diy::DiscreteBounds            Bounds;
// typedef diy::RegularGridLink           RGLink;
// typedef diy::RegularDecomposer<Bounds> Decomposer;

// one point
struct Pt
{
    diy::Point<double, 3>    coords;
};

// // whether a point is inside given bounds
// bool inside(const Pt& pt, const Bounds bounds)
// {
//     for (int i = 0; i < 3; i++)
//         if (pt.coords[i] < bounds.min[i] || pt.coords[i] >= bounds.max[i] - 1)
//             return false;
//     return true;
// }

// one end point of a particle trace segment
struct EndPt
{
    int  pid;                                // particle ID
    Pt   pt;                                 // end pointof the trace
    Pt   pt_hold;                            // holding place for p{xyz} during balancing
    std::vector<std::vector<double>> cell_data;              // for storing zMid, zTop, vertVelTop
    int  sid;                                // segment ID of this part of the trace
    int  nsteps;                             // number of steps this particle went so far
    double zLevelParticle;
    int glCellIdx;                          // global cell index for current point
     int predonly;                          // if predonly>0 : pt should be filtered out after load balancing; predonly=1 => load pt; predonly=2 => cellid pt
    int source_gid;                         // only set for the cellid pts for reconstructing links

    const double& operator [](int i) const;
    double& operator [](int i);

    EndPt();
       
    EndPt(struct Segment& s);                // extract the end point of a segment
};

// one segment of a particle trace (trajectory)
struct Segment
{
    int        pid;                          // particle ID
    vector<Pt> pts;                          // points along trace
    int        sid;                          // segment ID of this part of the trace
    int        nsteps;                       // number of particle steps at start of segment

    Segment();
    Segment(EndPt& p);                        // construct a segment from one point
       
    // // whether end point is inside given bounds
    // bool inside(const int lb[3], const int ub[3]) const
    //     {
    //         for (int i = 0; i < 3; i++)
    //             if (pts.back().coords[i] < lb[i] || pts.back().coords[i] >= ub[i] - 1)
    //                 return false;
    //         return true;
    //     }
};

struct Halo
{   
    int src_gid;
    int dest_gid;

    std::vector<int> glCellIDs;
    std::vector<double> xCell;
    std::vector<double> yCell;
    std::vector<double> zCell;

    std::vector<int> glVertexIDs;
    std::vector<double> xVertex, yVertex, zVertex;
    std::vector<double> velocityXv;
    std::vector<double> velocityYv;
    std::vector<double> velocityZv;
    std::vector<double> vertVelocityTop;
    std::vector<double> zMid;
    std::vector<double> zTop;

};



// specialize the serialization of a segment and Halo
namespace diy
{   

    template<>
    struct Serialization<EndPt>
    {
        static
        void save(diy::BinaryBuffer& bb, const EndPt& x)
            {
                       diy::save(bb, x.pid);
                       diy::save(bb, x.pt);
                       diy::save(bb, x.pt_hold);
                       diy::save(bb, x.cell_data);
                       diy::save(bb, x.sid);
                       diy::save(bb, x.nsteps);
                       diy::save(bb, x.zLevelParticle);
                       diy::save(bb, x.glCellIdx);
                       diy::save(bb, x.predonly);
                       diy::save(bb, x.source_gid);
            }
        static
        void load(diy::BinaryBuffer& bb, EndPt& x)
            {
                       diy::load(bb, x.pid);
                       diy::load(bb, x.pt);
                       diy::load(bb, x.pt_hold);
                       diy::load(bb, x.cell_data);
                       diy::load(bb, x.sid);
                       diy::load(bb, x.nsteps);
                       diy::load(bb, x.zLevelParticle);
                       diy::load(bb, x.glCellIdx);
                       diy::load(bb, x.predonly);
                       diy::load(bb, x.source_gid);
            }
    };

    template<>
    struct Serialization<Segment>
    {
        static
        void save(diy::BinaryBuffer& bb, const Segment& x)
            {
                diy::Serialization<int>::save(bb, x.pid);
                diy::Serialization< vector <Pt> >::
                    save(bb, static_cast< const vector<Pt>& >(x.pts));
                diy::Serialization<int>::save(bb, x.sid);
            }
        static
        void load(diy::BinaryBuffer& bb, Segment& x)
            {
                diy::Serialization<int>::load(bb, x.pid);
                diy::Serialization< vector<Pt> >::
                    load(bb, static_cast< vector<Pt>& >(x.pts));
                diy::Serialization<int>::load(bb, x.sid);
            }
    };

    template<>
     struct Serialization<Halo>
    {
        static
        void save(diy::BinaryBuffer& bb, const Halo& x)
            {
                diy::save(bb, x.src_gid);
                diy::save(bb, x.dest_gid);
                diy::save(bb, x.glCellIDs);
                diy::save(bb, x.xCell);
                diy::save(bb, x.yCell);
                diy::save(bb, x.zCell);
                diy::save(bb, x.glVertexIDs);
                diy::save(bb, x.velocityXv);
                diy::save(bb, x.velocityYv);
                diy::save(bb, x.velocityZv);
                diy::save(bb, x.vertVelocityTop);
                diy::save(bb, x.xVertex);
                diy::save(bb, x.yVertex);
                diy::save(bb, x.zVertex);
                diy::save(bb, x.zMid);
                diy::save(bb, x.zTop);
            }
        static
        void load(diy::BinaryBuffer& bb, Halo& x)
            {
                diy::load(bb, x.src_gid);
                diy::load(bb, x.dest_gid);
                diy::load(bb, x.glCellIDs);
                diy::load(bb, x.xCell);
                diy::load(bb, x.yCell);
                diy::load(bb, x.zCell);
                diy::load(bb, x.glVertexIDs);
                diy::load(bb, x.velocityXv);
                diy::load(bb, x.velocityYv);
                diy::load(bb, x.velocityZv);
                diy::load(bb, x.vertVelocityTop);
                diy::load(bb, x.xVertex);
                diy::load(bb, x.yVertex);
                diy::load(bb, x.zVertex);
                diy::load(bb, x.zMid);
                diy::load(bb, x.zTop);

            }
    };


}

#endif