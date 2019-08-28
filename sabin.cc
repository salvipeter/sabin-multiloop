// Based on:
//
// M. A. Sabin: Further transfinite developments
// In: The Mathematics of Surfaces VIII (Ed. R. Cripps), pp. 161-173, 1998.

#include <fstream>
#include <iostream>

#include <Eigen/Dense>

#include <geometry.hh>

using namespace Geometry;

// Global parameters
const double multiplier = 10.0;
const size_t resolution = 30;
const double domain_tolerance = 1.0e-5;

using Boundary =
  double(*)(double &u, double &v, Point3D &point, Vector3D &cross_derivative);


// Setup boundaries as in the paper (or at least similarly)

double b1(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  double t = u, d = v;
  v = 0;
  point = Point3D(0, 0, 0) * (1 - t) + Point3D(36, 0, 0) * t;
  cross_derivative = Vector3D(0, 1, 0) * multiplier;
  return d;
}

double b3(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  double t = u, d = 1 - v;
  v = 1;
  point = Point3D(0, 0, 20) * (1 - t) + Point3D(36, 0, 20) * t;
  cross_derivative = Vector3D(0, 1, 0) * multiplier;
  return d;
}

double b2(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  double phi = (1 - v) * M_PI, d = u;
  u = 0;
  point = Point3D(0, 0, 10) + Vector3D(0, sin(phi), cos(phi)) * 10;
  cross_derivative = Vector3D(1, 0, 0) * multiplier;
  return d;
}

double b4(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  double phi = (1 - v) * M_PI, d = 1 - u;
  u = 1;
  point = Point3D(36, 0, 10) + Vector3D(0, sin(phi), cos(phi)) * 10;
  cross_derivative = Vector3D(-1, 0, 0) * multiplier;
  return d;
}

double hole1(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  Point2D center(0.3, 0.5), p(u, v);
  double radius = 0.1;
  auto q = center + (p - center).normalize() * radius;
  u = q[0];
  v = q[1];
  point = Point3D(36 * q[0], 15, 20 * q[1]);
  cross_derivative = Vector3D(0, -1, 0) * multiplier;
  return (p - q).norm();
}

double hole2(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  Point2D center(0.7, 0.5), p(u, v);
  double radius = 0.1;
  auto q = center + (p - center).normalize() * radius;
  u = q[0];
  v = q[1];
  point = Point3D(36 * q[0], 15, 20 * q[1]);
  cross_derivative = Vector3D(0, -1, 0) * multiplier;
  return (p - q).norm();
}


// Main code

std::pair<Point3D, Vector3D>
evalPatch(const std::vector<Boundary> &boundaries, double u, double v) {
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero(), b = Eigen::Matrix3d::Zero();
  double ui, vi, di;
  Point3D Pi;
  Vector3D Ti;

  for (const auto &boundary : boundaries) {
    ui = u; vi = v;
    di = boundary(ui, vi, Pi, Ti);
    if (di < domain_tolerance)
      return { Pi, Vector3D(0, 0, 0) }; // kutykurutty

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
      b(1, k) += (3 * (u - ui) * Pi[k] + (u - ui) * di * Ti[k]) * denom;
      b(2, k) += (3 * (v - vi) * Pi[k] + (v - vi) * di * Ti[k]) * denom;
    }
  }

  Eigen::Matrix3d x = A.colPivHouseholderQr().solve(b);
  Point3D p(x(0,0), x(0,1), x(0,2));
  Vector3D j1(-x(1,0), -x(1,1), -x(1,2)), j2(-x(2,0), -x(2,1), -x(2,2));
  return { p, (j1 ^ j2).normalize() };
}

int main() {
  std::vector<Boundary> boundaries = { b1, b2, b3, b4 };
  std::ofstream f("/tmp/sabin.obj");
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    for (size_t j = 0; j <= resolution; ++j) {
      double v = (double)j / resolution;
      auto p = evalPatch(boundaries, u, v).first;
      f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    }
  }
  for (size_t i = 0; i < resolution; ++i)
    for (size_t j = 0; j < resolution; ++j) {
      size_t index = i * (resolution + 1) + j;
      f << "f " << index + 1
        << ' '  << index + 2
        << ' '  << index + resolution + 3 << std::endl;
      f << "f " << index + 1
        << ' '  << index + resolution + 3
        << ' '  << index + resolution + 2 << std::endl;
    }
}
