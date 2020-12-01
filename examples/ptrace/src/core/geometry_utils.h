#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H


#include <math.h>
#include <vector>
#include "misc.h"

class mpas_io;
class block;

double triangle_signed_area(double *a, double *b, double *c, double radius, int flag=0);

inline double get_arc_length(double ax, double ay, double az,
					double bx, double by, double bz, int flag=0){


	double cx = bx - ax;
    double cy = by - ay;
    double cz = bz - az;

    double r = sqrt(ax*ax + ay*ay + az*az);
    double c = sqrt(cx*cx + cy*cy + cz*cz);

    // if (flag==1){
    //     double b = sqrt(bx*bx + by*by + bz*bz);
    //     dprint("insine %f, r %f %f", c/(2.0*r), r, b);
    // }
    double mpas_arc_length = r * 2.0 * asin(c/(2.0*r));

	return mpas_arc_length;
}


int get_vertical_id(int nLevels, double zLoc, double *zMid);

void get_nearby_cell_index(int nCells,
                           const double *xc, 
                           const double *yc, 
                           const double *zc,
                           const double xp, 
                           const double yp, 
                           const double zp,
                           const block &mpas1, 
                           int &lastCell, 
                           const int *cellsOnCell,
                           const int *nEdgesOnCell
                           );

void get_nearby_cell_index_sl(int nCells,
                           const double *xc, 
                           const double *yc, 
                           const double *zc,
                           const double xp, 
                           const double yp, 
                           const double zp,
                           const block &mpas1, 
                           int &lastCell, 
                           const int *cellsOnCell,
                           const int *nEdgesOnCell
                           );


#endif // GEOMETRY_UTILS