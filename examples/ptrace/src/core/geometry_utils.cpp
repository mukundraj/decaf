#include "geometry_utils.h"
#include <math.h>
#include <algorithm> 
#include "misc.h"
#include "block.h"
#include <Eigen/Dense>

double triangle_signed_area(double *a, double *b, double *c, double radius){

  double ab = get_arc_length(a[0], a[1], a[2], b[0], b[1], b[2])/radius;
  double bc = get_arc_length(b[0], b[1], b[2], c[0], c[1], c[2])/radius;
  double ca = get_arc_length(c[0], c[1], c[2], a[0], a[1], a[2])/radius;
  double semiperim = 0.5 * (ab + bc + ca);


  double tanqe = sqrt(std::max(0.0,tan(0.5* semiperim) * tan(0.5* (semiperim - ab))
   * tan(0.5 * (semiperim - bc)) * tan(0.5 * (semiperim - ca))));


 double mpas_triangle_signed_area_sphere = 4.0 * radius * radius * atan(tanqe);

      // ! computing correct signs (in similar fashion to mpas_sphere_angle)

 double ablen[3], aclen[3], dlen[3];
 ablen[1] = b[1] - a[1];
 ablen[2] = b[2] - a[2];
 ablen[0] = b[0] - a[0];

 aclen[1] = c[1] - a[1];
 aclen[2] = c[2] - a[2];
 aclen[0] = c[0] - a[0];

 dlen[0] =   (ablen[1] * aclen[2]) - (ablen[2] * aclen[1]);
 dlen[1] = -((ablen[0] * aclen[2]) - (ablen[2] * aclen[0]));
 dlen[2] =   (ablen[0] * aclen[1]) - (ablen[1] * aclen[0]);

 if ((dlen[0]*a[0] + dlen[1]*a[1] + dlen[2]*a[2]) < 0.0) 
  mpas_triangle_signed_area_sphere = -mpas_triangle_signed_area_sphere;


return mpas_triangle_signed_area_sphere;
}

int get_vertical_id(int nLevels, double zLoc, double *zMid){

  // dprint("nLevels %d", nLevels);
  for (int aLevel=0; aLevel<nLevels-1; aLevel++){

    if (zMid[aLevel+1]<=zLoc && zLoc <= zMid[aLevel]){
      // the point is bounded by the levels (store the bottom level)
      return aLevel;
    }

  }


  if (zLoc < zMid[nLevels-1]){
    // case where location is smallest value
    return -1;
  }
  else if (zLoc > zMid[0]){
    // case where location is largest value
    return 0;
  }

}

// return global cell id with fortran starting index
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
                           const int *nEdgesOnCell)
{

  double pointRadius, xPoint[3];
  std::map<int, Eigen::Array3d> xCell;
  int aPoint;
  int cellID;

 
  
    int cellGuess = lastCell;
    cellID = -1;

    while (cellID != cellGuess)
    {
      
      // we have a known cell
      cellID = cellGuess;

      // normalize locations to same spherical shell (unit) for direct comparison

      // for point itself
      pointRadius = sqrt(xp * xp + yp * yp + zp * zp);
      xPoint[0] = xp / pointRadius;
      xPoint[1] = yp / pointRadius;
      xPoint[2] = zp / pointRadius;

      int localCellID = cellID;
      xCell[cellGuess] <<  xc[localCellID]/pointRadius, yc[localCellID]/pointRadius, zc[localCellID]/pointRadius;
      
      // for point neighbors
      int aPoint_local = -7;

      for (int iPoint = 0; iPoint < nEdgesOnCell[localCellID]; iPoint++)
      {
        
        aPoint = cellsOnCell[localCellID * mpas1.maxEdges + iPoint]; // global cell id of neighbor
        aPoint = aPoint - 1;
        // aPoint_ = aPoint;
        if (aPoint > nCells || aPoint == 0)
          continue; // todo: confirm logic

        // aPoint_local = aPoint - 1; // local cell id of neighbor
        aPoint_local = aPoint; // local cell id of neighbor
        pointRadius = mpas1.cradius;
        xCell[aPoint] <<  xc[aPoint_local]/pointRadius, yc[aPoint_local]/pointRadius, zc[aPoint_local]/pointRadius;
       
      }

        double dx = xPoint[0] - xCell[cellID](0);
        double dy = xPoint[1] - xCell[cellID](1);
        double dz = xPoint[2] - xCell[cellID](2);
        double r2Min = dx*dx + dy*dy + dz*dz;
        double r2;

        // dprint("r2Min %.12e", r2Min);
        
        for (int iPoint=0; iPoint<nEdgesOnCell[localCellID]; ++iPoint){
          aPoint = cellsOnCell[localCellID*mpas1.maxEdges +iPoint] -1 ;

          // if (aPoint > nCells || aPoint == 0)
          if (aPoint > nCells)
            continue; // todo: confirm logic

          // compute squared distances

          dx = xPoint[0] - xCell[aPoint](0);
          dy = xPoint[1] - xCell[aPoint](1);
          dz = xPoint[2] - xCell[aPoint](2);
          r2 = dx*dx + dy*dy + dz*dz;
          // dprint("r2 %.12e", r2);
          if ( r2 < r2Min){
            // we have a new closest point
            cellGuess = aPoint;
            r2Min = r2;
          }
        }
      
    }

      lastCell = cellID;
}


// return global cell id with fortran starting index
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
                           const int *nEdgesOnCell)
{

  double pointRadius, xPoint[3];
  std::map<int, Eigen::Array3d> xCell;
  int aPoint;
  int cellID;

  if (lastCell < 1)
  {
    // if (1){
    // brute force solution
    Eigen::VectorXd q(3);
    q << xp, yp, zp;
    // get nearest cell neighbor ids using mpas_c
    Eigen::VectorXi nearest_cell_idx(1);
    Eigen::VectorXd dists2_cell(1);
    mpas1.nns_cells->knn(q, nearest_cell_idx, dists2_cell, 1, 0, Nabo::NNSearchF::SORT_RESULTS | Nabo::NNSearchF::ALLOW_SELF_MATCH);
    cellID = nearest_cell_idx[0] + 1;
    // dprint("cellID %d", cellID);
  }
  else
  {

    int cellGuess = lastCell;
    cellID = -1;

    while (cellID != cellGuess)
    {

      // we have a known cell
      cellID = cellGuess;

      // normalize locations to same spherical shell (unit) for direct comparison

      // for point itself
      pointRadius = sqrt(xp * xp + yp * yp + zp * zp);
      xPoint[0] = xp / pointRadius;
      xPoint[1] = yp / pointRadius;
      xPoint[2] = zp / pointRadius;

      // dprint("cellID %d", cellID);
      int localCellID = mpas1.cellIndex.at(cellID);
      xCell[cellGuess] <<  xc[localCellID]/pointRadius, yc[localCellID]/pointRadius, zc[localCellID]/pointRadius;
      
      // for point neighbors
      int aPoint_local = -7;

      for (int iPoint = 0; iPoint < nEdgesOnCell[localCellID]; iPoint++)
      {

        // aPoint = cellsOnCell[cellID*mpas1.maxEdges+iPoint] - 1; // local cell id of neighbor
        // iPoint_ = iPoint;
        // aPoint = cellsOnCell[localCellID * mpas1.maxEdges + iPoint] - 1; // local cell id of neighbor
        aPoint = cellsOnCell[localCellID * mpas1.maxEdges + iPoint]; // global cell id of neighbor

        // dprint("aPoint %d", aPoint);
        
       
        // aPoint_ = aPoint;
        if (aPoint > nCells || aPoint == 0)
          continue; // todo: confirm logic

        // dprint("aPoint %d", aPoint);
        aPoint_local = mpas1.cellIndex.at(aPoint); // local cell id of neighbor
        pointRadius = mpas1.cradius;
        xCell[aPoint] <<  xc[aPoint_local]/pointRadius, yc[aPoint_local]/pointRadius, zc[aPoint_local]/pointRadius;
        // pointRadius = sqrt(xc[aPoint_local] * xc[aPoint_local] + yc[aPoint_local] * yc[aPoint_local] + zc[aPoint_local] * zc[aPoint_local]);

        // dprint("xCell[%d], %f %f %f, locCellID %d xc %f pointR %f", aPoint, xCell[aPoint](0), xCell
        // [aPoint](1), xCell[aPoint](2), localCellID, xc[localCellID], pointRadius);
      }
      // dprint("--");

      // dprint("cellID %d cellGuess %d", cellID, cellGuess);

      
        double dx = xPoint[0] - xCell[cellID](0);
        double dy = xPoint[1] - xCell[cellID](1);
        double dz = xPoint[2] - xCell[cellID](2);
        double r2Min = dx*dx + dy*dy + dz*dz;
        double r2;

        // dprint("r2Min %.12e", r2Min);
        
        for (int iPoint=0; iPoint<nEdgesOnCell[localCellID]; ++iPoint){
          aPoint = cellsOnCell[localCellID*mpas1.maxEdges+iPoint];

          if (aPoint > nCells || aPoint == 0)
            continue; // todo: confirm logic
          aPoint -= 1;
          // compute squared distances

          dx = xPoint[0] - xCell[aPoint](0);
          dy = xPoint[1] - xCell[aPoint](1);
          dz = xPoint[2] - xCell[aPoint](2);
          r2 = dx*dx + dy*dy + dz*dz;
          // dprint("r2 %.12e", r2);
          if ( r2 < r2Min){
            // we have a new closest point
            cellGuess = aPoint;
            r2Min = r2;
          }
        }
        
    }
  }

  lastCell = cellID;
}
