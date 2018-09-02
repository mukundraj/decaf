#ifndef PATHLINES_H
#define PATHLINES_H

#include <vector>
#include <string>
#include "flow.h"
#include <Eigen/Dense>
#include "mpaso.h"
#include "nabo/nabo.h"


class pathlines
{



 public:

  pathlines(std::string &ip_file_name, std::string &op_file_name,
              unsigned int init_t, unsigned int fin_t, unsigned int interval, float h, std::vector<double> &seeds_xyz);

  void runge_kutta(double &cx, double &cy, double &cz, float h, Eigen::VectorXi &nearest_idx, bool &terminate_stat);
  void get_curpos_velocity_pl(double cx, double cy, double cz, Eigen::Vector3d &c_vel, Eigen::VectorXi &nearest_idx, bool &terminate_stat);
  void generate_seeds(mpaso &mpasobj, int skipval, std::vector<double> &seeds_xyz);

  ~pathlines();

  void compute();
  void write();


  mpaso mpas1, mpas2;
//  float cur_t, fir_t, sec_t;
  Eigen::Vector3d c_vel={0,0,0}; // placeholder for velocity values
  const double radius = 6371220.;

  Nabo::NNSearchD*  nns_cells;
  Nabo::NNSearchD*  nns;

  unsigned int fir_t, sec_t, cur_t;


};



#endif
