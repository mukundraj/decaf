#ifndef INTERPOLATORS_H
#define INTERPOLATORS_H

#include <cstdlib>
#include <vector>
#include <Eigen/Dense>

template <typename T>
static inline void cross_product(const T A[3], const T B[3], T C[3])
{
  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];
}

template <typename T>
static inline T dot_product(const T A[3], const T B[3])
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template <typename T>
static inline T barycentric_point2triangle(const T P0[3], const T P1[3], const T P2[3], const T P[3], T lambda[3])
{
  T u[3] = {P1[0] - P0[0], P1[1] - P0[1], P1[2] - P0[2]},
    v[3] = {P2[0] - P0[0], P2[1] - P0[1], P2[2] - P0[2]},
    w[3] = {P[0] - P0[0], P[1] - P0[1], P[2] - P0[2]};
  T n[3], uw[3], wv[3];

  cross_product(u, v, n); // n = u x v
  cross_product(u, w, uw); // uw = u x w
  cross_product(w, v, wv); // wv = w x v

  T n2 = dot_product(n, n); // n^2

  lambda[2] = dot_product(uw, n) / n2;
  lambda[1] = dot_product(wv, n) / n2;
  lambda[0] = 1 - lambda[2] - lambda[1];
}

void interpolate_vertically(size_t nVertLevels, std::vector<double> &zTopVertex, Eigen::VectorXi &nearest_idx,
                            std::vector<double> &values, double depth,
                            std::vector<double> &xVertex, std::vector<double> &yVertex, std::vector<double> &zVertex,
                            std::vector<double> &velocityXv, std::vector<double> &velocityYv, std::vector<double> &velocityZv, int gid);


double linear_inter(double x, double x1, double x2, double q00, double q01);

bool interpolate_horizontally(double cx, double cy, double cz,
                              std::vector<double> &values, Eigen::Vector3d &c_vel);

int bin_search_index(std::vector<double> &zTopVertex, int lo, int hi, double val);

#endif
