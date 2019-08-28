#include "solver.hh"

#include <Eigen/Dense>

const double domain_tolerance = 1.0e-5;

std::pair<Point3D, Vector3D>
evalPatch(const std::vector<Boundary> &boundaries, double u, double v) {
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero(), b = Eigen::Matrix3d::Zero();
  double ui, vi, di;
  Point3D Pi;
  Vector3D Ti, du, dv;

  for (const auto &boundary : boundaries) {
    ui = u; vi = v;
    di = boundary(ui, vi, Pi, du, dv);
    if (di < domain_tolerance)
      return { Pi, (du ^ dv).normalize() };
    Ti = du * (u - ui) + dv * (v - vi); // = cross derivative scaled by di

    double denom = std::pow(di, -3);
    A(0, 0) += 2 * denom;
    A(0, 1) += (u - ui) * denom;
    A(0, 2) += (v - vi) * denom;
    A(1, 0) += 3 * (u - ui) * denom;
    A(1, 1) += 2 * std::pow(u - ui, 2) * denom;
    A(1, 2) += 2 * (u - ui) * (v - vi) * denom;
    A(2, 0) += 3 * (v - vi) * denom;
    A(2, 1) += 2 * (u - ui) * (v - vi) * denom;
    A(2, 2) += 2 * std::pow(v - vi, 2) * denom;

    for (size_t k = 0; k < 3; ++k) {
      b(0, k) += (2 * Pi[k] + Ti[k]) * denom;
      b(1, k) += (3 * (u - ui) * Pi[k] + (u - ui) * Ti[k]) * denom;
      b(2, k) += (3 * (v - vi) * Pi[k] + (v - vi) * Ti[k]) * denom;
    }
  }

  Eigen::Matrix3d x = A.colPivHouseholderQr().solve(b);
  Point3D p(x(0,0), x(0,1), x(0,2));
  Vector3D j1(-x(1,0), -x(1,1), -x(1,2)), j2(-x(2,0), -x(2,1), -x(2,2));
  return { p, (j1 ^ j2).normalize() };
}
