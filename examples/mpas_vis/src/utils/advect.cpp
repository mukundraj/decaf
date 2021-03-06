//---------------------------------------------------------------------------
//
// particle advection
//
// courtesy Hanqi Guo
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
//
// copied with permission by Tom Peterka
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------

#include "advect.h"
#include "lerp.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>

// extern "C" {

    bool trace_3D_brown(const int *st,
                        const int *sz,
                        const float **vec,
                        float *X,
                        float h,
                        float *Y = NULL)     // result returned here if not NULL, otherwise in X
    {
        if (!inside(3, st, sz, X)) return false;

        if (Y)
        {
            Y[0] = X[0] + (float)rand()/RAND_MAX - 0.5;
            Y[1] = X[1] + (float)rand()/RAND_MAX - 0.5;
            Y[2] = X[2] + (float)rand()/RAND_MAX - 0.5;
        }
        else
        {
            X[0] = X[0] + (float)rand()/RAND_MAX - 0.5;
            X[1] = X[1] + (float)rand()/RAND_MAX - 0.5;
            X[2] = X[2] + (float)rand()/RAND_MAX - 0.5;
        }

        return true;
    }

   //  bool trace_3d_rk1_mpas(mpaso &mpas_g,
   //                                  mpaso &mpas_c,
   //                              double *X,
   //                              double h,
   //                              double *Y)
   // {
   //
   //          return false;
   //
   //  }

    bool trace_3D_rk1(const int *gst,
                      const int *gsz,
                      const int *st,
                      const int *sz,
                      const float **vec,
                      float *X,
                      float h,
                      float *Y = NULL)       // result returned here if not NULL, otherwise in X
    {
        if (!inside(3, gst, gsz, X)) return false;

        float v[3];
        if (!lerp3D(X, gst, gsz, 3, vec, v))
            return false;

        if (Y)
        {
            Y[0] = X[0] + h*v[0];
            Y[1] = X[1] + h*v[1];
            Y[2] = X[2] + h*v[2];
        }
        else
        {
            X[0] = X[0] + h*v[0];
            X[1] = X[1] + h*v[1];
            X[2] = X[2] + h*v[2];
        }

        return true;
    }

    bool trace_4D_rk1(const int *gst,
                      const int *gsz,
                      const int *st,
                      const int *sz,
                      const float **vec,
                      float *X,
                      float h,
                      float *Y = NULL)       // result returned here if not NULL, otherwise in X
    {
        if (!inside(4, gst, gsz, X)) return false;

        float v[3];
        if (!lerp4D(X, gst, gsz, 3, vec, v))
            return false;

        if (Y)
        {
            Y[0] = X[0] + h*v[0];
            Y[1] = X[1] + h*v[1];
            Y[2] = X[2] + h*v[2];
            Y[3] = X[3] + h; // TODO
        }
        else
        {
            X[0] = X[0] + h*v[0];
            X[1] = X[1] + h*v[1];
            X[2] = X[2] + h*v[2];
            X[3] = X[3] + h; // TODO
        }

        return true;
    }

// }
