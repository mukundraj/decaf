//---------------------------------------------------------------------------
//
// tessellation
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// mdreher@anl.gov
//
//--------------------------------------------------------------------------

#include <decaf/decaf.hpp>
#include <decaf/data_model/simplefield.hpp>
#include <decaf/data_model/arrayfield.hpp>
#include <decaf/data_model/blockfield.hpp>
#include <decaf/data_model/boost_macros.h>

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <map>
#include <cstdlib>

#include "wflow.hpp"                         // defines the workflow for this example
#include "tess/tess.h"
#include "tess/tess.hpp"

// using namespace decaf;
using namespace std;

// consumer
void tessellate(Decaf* decaf, MPI_Comm comm)
{
    // this first test just creates and tessellates some synthetic data to demostrate that
    // tess (a diy program) can be run as a decaf node task

    int tot_blocks  = 8;                      // total number of blocks in the domain
    int mem_blocks  = -1;                     // max blocks in memory
    int dsize[3]    = {11, 11, 11};           // domain grid size
    float jitter    = 2.0;                    // max amount to randomly displace particles
    int wrap        = 0;                      // wraparound neighbors flag
    int num_threads = 1;                      // threads diy can use
    char outfile[256];                        // output file name
    strcpy(outfile, "del.out");
    double times[TESS_MAX_TIMES];             // timing

    // data extents
    typedef     diy::ContinuousBounds         Bounds;
    Bounds domain;
    for(int i = 0; i < 3; i++)
    {
        domain.min[i] = 0;
        domain.max[i] = dsize[i] - 1.0;
    }

    // event loop
    // TODO: verify that all of the below can run iteratively, only tested for one time step
    vector<pConstructData> in_data;
    while (decaf->get(in_data))
    {

        ArrayFieldf xyzpos = in_data[0]->getFieldData<ArrayFieldf>("xyz_pos");
        if(!xyzpos)
        {
            fprintf(stderr, "ERROR tessellate: unable to find the field required \"xyz_pos\" "
                    "in the data model.\n");
            return;
        }

        BlockField block = in_data[0]->getFieldData<BlockField>("domain_block");
        if(!block)
        {
            fprintf(stderr, "ERROR tessellate: unable to find the field required \"domain_block\" "
                    "in the data model.\n");
            return;
        }

        float* xyz = xyzpos.getArray();
        Block<3>* blk = block.getBlock();
        float* global_box = blk->getGlobalBBox(); // min and size (not min and max)
        float* local_box = blk->getLocalBBox();   // min and size (not min and size)

        // TODO: use the positions I got from decaf
        // don't free anything that I get from decaf, its smart pointers take care of it

        // debug
        int nbParticle = xyzpos->getNbItems();
        // fprintf(stderr, "nbParticle=%d\n", nbParticle);
        // for (int i = 0; i < nbParticle; i++)
        //     fprintf(stderr, "%.3f %.3f %.3f\n", xyz[3*i], xyz[3*i+1], xyz[3*i+2]);

        // init diy
        diy::mpi::communicator    world(comm);
        diy::FileStorage          storage("./DIY.XXXXXX");
        diy::Master               master(world,
                                         num_threads,
                                         mem_blocks,
                                         &create_block,
                                         &destroy_block,
                                         &storage,
                                         &save_block,
                                         &load_block);
        diy::RoundRobinAssigner   assigner(world.size(), tot_blocks);
        AddAndGenerate            create(master, jitter);

        // decompose
        std::vector<int> my_gids;
        assigner.local_gids(world.rank(), my_gids);
        diy::RegularDecomposer<Bounds>::BoolVector          wraps;
        diy::RegularDecomposer<Bounds>::BoolVector          share_face;
        diy::RegularDecomposer<Bounds>::CoordinateVector    ghosts;
        if (wrap)
            wraps.assign(3, true);
        diy::decompose(3, world.rank(), domain, assigner, create, share_face, wraps, ghosts);

        // tessellate
        quants_t quants;
        timing(times, -1, -1, world);
        timing(times, TOT_TIME, -1, world);
        tess(master, quants, times);

        // output
        tess_save(master, outfile, times);
        timing(times, -1, TOT_TIME, world);
        tess_stats(master, quants, times);

        // send a token downstream that data are ready
        // TODO: send data instead of the token (don't write to a file)
        fprintf(stderr, "tess sending token\n");
        int value = 0;

        SimpleFieldi data(value);
        pConstructData container;
        container->appendData(string("var"), data,
                              DECAF_NOFLAG, DECAF_PRIVATE,
                              DECAF_SPLIT_KEEP_VALUE, DECAF_MERGE_ADD_VALUE);
        decaf->put(container);
    } // decaf event loop

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    fprintf(stderr, "tessellation terminating\n");
    decaf->terminate();
}

int main(int argc,
         char** argv)
{
    // define the workflow
    Workflow workflow;
    make_wflow(workflow);

    MPI_Init(NULL, NULL);

    // create decaf
    Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

    // start the task
    tessellate(decaf, decaf->con_comm_handle());

    // cleanup
    delete decaf;
    MPI_Finalize();

    fprintf(stderr, "finished prod\n");
    return 0;
}
