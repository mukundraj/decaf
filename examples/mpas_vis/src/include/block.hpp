//---------------------------------------------------------------------------
//
// diy2-vtk7 parallel particle advection block class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------
#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi/datatypes.hpp>
#include <diy/io/bov.hpp>
#include <diy/pick.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>

#ifdef WITH_VTK
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#endif

#include <fstream>
#include <stdio.h>

//typedef diy::DiscreteBounds            Bounds;
//typedef diy::RegularGridLink           RGLink;
//typedef diy::RegularDecomposer<Bounds> Decomposer;

using namespace std;

// the diy block
struct DBlock
{
    DBlock() : nvecs(0), init(0), done(0) {}
    ~DBlock()
    {
        if (nvecs)
        {
            for (int i = 0; i < 3; i++)
                delete vel[i];
        }
    }

    static void* create()
    {
        return new DBlock;
    }
    static void destroy (void* b)
    {
        delete static_cast<DBlock*>(b);
    }
    static void save(const void* b_, diy::BinaryBuffer& bb)
    {
        const DBlock* b = static_cast<const DBlock*>(b_);
        diy::save(bb, b->nvecs);
        diy::save(bb, b->vel[0], b->nvecs);
        diy::save(bb, b->vel[1], b->nvecs);
        diy::save(bb, b->vel[2], b->nvecs);
        diy::save(bb, b->init);
        diy::save(bb, b->done);
        // TODO: serialize vtk structures
    }
    static void load(void* b_, diy::BinaryBuffer& bb)
    {
        DBlock* b = static_cast<DBlock*>(b_);
        diy::load(bb, b->nvecs);
        b->vel[0] = new float[b->nvecs];
        b->vel[1] = new float[b->nvecs];
        b->vel[2] = new float[b->nvecs];
        diy::load(bb, b->vel[0], b->nvecs);
        diy::load(bb, b->vel[1], b->nvecs);
        diy::load(bb, b->vel[2], b->nvecs);
        diy::load(bb, b->init);
        diy::load(bb, b->done);
        // TODO: serialize vtk structures
    }


    float                *vel[3];            // pointers to vx, vy, vz arrays (v[0], v[1], v[2])
    size_t               nvecs;              // number of velocity vectors
    int                  init, done;         // initial and done flags
    vector<Segment> segments;                // finished segments of particle traces
};
