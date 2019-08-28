// Based on:
//
// M. A. Sabin: Further transfinite developments
// In: The Mathematics of Surfaces VIII (Ed. R. Cripps), pp. 161-173, 1998.

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>

#include "solver.hh"

// Global parameters
const double multiplier = 30.0;
const size_t resolution = 100;


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
  bool inside = (p - center).norm() < radius;
  auto q = center + (p - center).normalize() * radius;
  u = q[0];
  v = q[1];
  point = Point3D(36 * q[0], 15, 36 * q[1] - 8);
  cross_derivative = Vector3D(0, -1, 0) * multiplier;
  return inside ? -1 : (p - q).norm();
}

double hole2(double &u, double &v, Point3D &point, Vector3D &cross_derivative) {
  Point2D center(0.7, 0.5), p(u, v);
  double radius = 0.1;
  bool inside = (p - center).norm() < radius;
  auto q = center + (p - center).normalize() * radius;
  u = q[0];
  v = q[1];
  point = Point3D(36 * q[0], 15, 36 * q[1] - 8);
  cross_derivative = Vector3D(0, -1, 0) * multiplier;
  return inside ? -1 : (p - q).norm();
}


// Main code

int main() {
  std::vector<Boundary> boundaries = { b1, b2, b3, b4 }, holes = { hole1, hole2 };

  std::set<size_t> hole_indices;
  auto is_in_hole = [&](size_t i) { return hole_indices.find(i) != hole_indices.end(); };
  boundaries.insert(boundaries.end(), holes.begin(), holes.end());

  std::ofstream f("/tmp/sabin.obj");
  size_t index = 0;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    for (size_t j = 0; j <= resolution; ++j, ++index) {
      double v = (double)j / resolution;

      // Heuristic handling of holes
      bool in_hole = false;
      for (const auto &hole : { hole1, hole2 }) {
        double u1 = u, v1 = v;
        Point3D p;
        Vector3D t;
        if (hole(u1, v1, p, t) < 0) {
          in_hole = true;
          hole_indices.insert(index);
          f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
          break;
        }
      }
      if (in_hole)
        continue;

      auto p = evalPatch(boundaries, u, v).first;
      f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    }
  }
  for (size_t i = 0; i < resolution; ++i)
    for (size_t j = 0; j < resolution; ++j) {
      index = i * (resolution + 1) + j + 1;
      if (!is_in_hole(index) || 
          !is_in_hole(index + 1) ||
          !is_in_hole(index + resolution + 2))
        f << "f " << index
          << ' '  << index + 1
          << ' '  << index + resolution + 2 << std::endl;
      if (!is_in_hole(index) ||
          !is_in_hole(index + resolution + 2) ||
          !is_in_hole(index + resolution + 1))
        f << "f " << index
          << ' '  << index + resolution + 2
          << ' '  << index + resolution + 1 << std::endl;
    }
}
