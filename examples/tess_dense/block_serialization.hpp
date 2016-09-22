//---------------------------------------------------------------------------
//
// diy block serialization (tess data structure)
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// mdreher@anl.gov
//
//--------------------------------------------------------------------------
#ifndef _BLOCK_SERIALIZATION
#define _BLOCK_SERIALIZATION

#include <decaf/data_model/baseconstructdata.hpp>
#include <boost/serialization/array.hpp>
#include "tess/tess.h"
#include "tess/tess.hpp"

using namespace std;

// block to be serialized
// derived from tess dblock_t with some fields removed and others added
struct SerBlock {
    int           gid;                // global block id
    float         mins[3], maxs[3];   // block extents
    struct bb_c_t data_bounds;        // global data extents
    struct bb_c_t box;	              // current box used in swap-reduce() when distributing particles
    int           num_orig_particles; // number of original particles
    int           num_particles;      // current number of particles
    float*        particles;          // all particles, original plus those received from neighbors
    int           num_tets;           // number of delaunay tetrahedra
    struct tet_t* tets;               // delaunay tets
    int*          rem_gids;           // owners of remote particles
    int*          rem_lids;	      // local ids of the remote particles
    int*          vert_to_tet;        // a tet that contains the vertex
    int           complete;           // whether tessellation is done
    vector<diy::BlockID> neighbors;   // block neighbors
};

BOOST_SERIALIZATION_SPLIT_FREE(SerBlock)
BOOST_IS_BITWISE_SERIALIZABLE(bb_c_t)
BOOST_IS_BITWISE_SERIALIZABLE(tet_t)
BOOST_IS_BITWISE_SERIALIZABLE(gb_t)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        void serialize(Archive& ar, bb_c_t& b, const unsigned int)
        {
            ar & b.min;
            ar & b.max;
        }
        template<class Archive>
        void save(Archive& ar, const SerBlock& b, unsigned int)
        {
            ar & BOOST_SERIALIZATION_NVP(b.gid);
            ar & BOOST_SERIALIZATION_NVP(b.mins);
            ar & BOOST_SERIALIZATION_NVP(b.maxs);
            ar & BOOST_SERIALIZATION_NVP(b.box);
            ar & BOOST_SERIALIZATION_NVP(b.data_bounds);
            ar & BOOST_SERIALIZATION_NVP(b.num_orig_particles);
            ar & BOOST_SERIALIZATION_NVP(b.num_particles);
            ar & boost::serialization::make_array<float>(b.particles, 3 * b.num_particles);
            ar & boost::serialization::make_array<int>(b.rem_gids,
                                                       b.num_particles - b.num_orig_particles);
            // TODO: this crashes, don't know why
            // ar & boost::serialization::make_array<int>(b.rem_lids,
            //                                            b.num_particles - b.num_orig_particles);
            ar & BOOST_SERIALIZATION_NVP(b.complete);
            ar & BOOST_SERIALIZATION_NVP(b.num_tets);
            ar & boost::serialization::make_array<tet_t>(b.tets, b.num_tets);
            ar & boost::serialization::make_array<int>(b.vert_to_tet, b.num_particles);
            ar & BOOST_SERIALIZATION_NVP(b.neighbors);
        }
        template<class Archive>
        void load(Archive& ar, SerBlock& b, unsigned int)
        {
            ar & BOOST_SERIALIZATION_NVP(b.gid);
            ar & BOOST_SERIALIZATION_NVP(b.mins);
            ar & BOOST_SERIALIZATION_NVP(b.maxs);
            ar & BOOST_SERIALIZATION_NVP(b.box);
            ar & BOOST_SERIALIZATION_NVP(b.data_bounds);
            ar & BOOST_SERIALIZATION_NVP(b.num_orig_particles);
            ar & BOOST_SERIALIZATION_NVP(b.num_particles);
            b.particles = new float[3 * b.num_particles];
            ar & boost::serialization::make_array<float>(b.particles, 3 * b.num_particles);
            b.rem_gids = new int[b.num_particles - b.num_orig_particles];
            ar & boost::serialization::make_array<int>(b.rem_gids,
                                                       b.num_particles - b.num_orig_particles);
            b.rem_lids = new int[b.num_particles - b.num_orig_particles];
            // TODO: this crashes, don't know why
            // ar & boost::serialization::make_array<int>(b.rem_lids,
            //                                            b.num_particles - b.num_orig_particles);
            ar & BOOST_SERIALIZATION_NVP(b.complete);
            ar & BOOST_SERIALIZATION_NVP(b.num_tets);
            b.tets = new tet_t[b.num_tets];
            ar & boost::serialization::make_array<tet_t>(b.tets, b.num_tets);
            b.vert_to_tet = new int[b.num_particles];
            ar & boost::serialization::make_array<int>(b.vert_to_tet, b.num_particles);
            ar & BOOST_SERIALIZATION_NVP(b.neighbors);
        }
    }
}

BOOST_CLASS_EXPORT_GUID(decaf::ArrayConstructData<SerBlock>, "ArrayConstructData<SerBlock>")

#endif