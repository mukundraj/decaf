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

#ifndef _ADVECT_H
#define _ADVECT_H

#include <stdbool.h>
#include "mpaso.h"
#ifdef __cplusplus
// extern "C" {
#endif

    bool trace_3D_rk1(const int   *gst,
                      const int   *gsz,
                      const int   *st,
                      const int   *sz,
                      const float **vec,
                      float       *X,
                      float        h,
                      float       *Y);

    bool trace_4D_rk1(const int   *gst,
                      const int   *gsz,
                      const int   *st,
                      const int   *sz,
                      const float **vec,
                      float       *X,
                      float        h,
                      float       *Y);

    bool trace_3D_brown(const int   *st,
                        const int   *sz,
                        const float **vec,
                        float       *X,
                        float        h,
                        float       *Y);

    bool trace_3D_rk1_mpas(mpaso &mpas_g, 
    				mpaso &mpas_c, 
				double *X, 
				double h, 
				double *Y);

#ifdef __cplusplus
// }
#endif

#endif
