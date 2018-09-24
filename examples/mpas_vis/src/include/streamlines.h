#ifndef STREAMLINES_H
#define STREAMLINES_H

#include <vector>
#include <string>
#include "flow.h"
#include <Eigen/Dense>
#include "mpaso.h"
#include "nabo/nabo.h"

#include <diy/assigner.hpp>



struct PBlock;

class streamlines
{



 public:
  streamlines();

  streamlines(mpaso &mpas1, std::vector<double> decaf_bar, std::string &op_file_name,
              float init_t, float fin_t, float h, std::vector<double> &seeds_xyz);

  void runge_kutta(double &cx, double &cy, double &cz, float h, Eigen::VectorXi &nearest_idx, bool &terminate_stat);
  void get_curpos_velocity_sl(double cx, double cy, double cz, Eigen::Vector3d &c_vel, Eigen::VectorXi &nearest_idx, bool &terminate_stat);
  void generate_seeds(mpaso &mpasobj, int skipval, std::vector<double> &seeds_xyz);

  ~streamlines();

  void compute();
  void write();


  mpaso *mpas1;
  float cur_t, fir_t, sec_t;
  Eigen::Vector3d c_vel={0,0,0}; // placeholder for velocity values
  const double radius = 6371220.;

  Nabo::NNSearchD*  nns_cells;
  Nabo::NNSearchD*  nns;

};



#endif
