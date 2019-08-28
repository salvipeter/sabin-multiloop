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

double b1(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv) {
  double t = u, d = v;
  v = 0;
  point = Point3D(0, 0, 0) * (1 - t) + Point3D(33, 0, 0) * t;
  du = Vector3D(1, 0, 0) * multiplier;
  dv = Vector3D(0, 1, 0) * multiplier;
  return d;
}

double b3(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv) {
  double t = u, d = 1 - v;
  v = 1;
  point = Point3D(0, 0, 20) * (1 - t) + Point3D(33, 0, 20) * t;
  du = Vector3D(1, 0, 0) * multiplier;
  dv = Vector3D(0, -1, 0) * multiplier;
  return d;
}

double b2(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv) {
  double phi = (1 - v) * M_PI, d = u;
  u = 0;
  point = Point3D(0, 0, 10) + Vector3D(0, sin(phi), cos(phi)) * 10;
  du = Vector3D(1, 0, 0) * multiplier;
  dv = Vector3D(0, -cos(phi), sin(phi)) * multiplier;
  return d;
}

double b4(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv) {
  double phi = (1 - v) * M_PI, d = 1 - u;
  u = 1;
  point = Point3D(33, 0, 10) + Vector3D(0, sin(phi), cos(phi)) * 10;
  du = Vector3D(1, 0, 0) * multiplier;
  dv = Vector3D(0, -cos(phi), sin(phi)) * multiplier;
  return d;
}

double hole1(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv) {
  Point2D center(0.3, 0.5), p(u, v);
  double radius = 0.1;
  auto deviation = p - center;
  bool inside = deviation.norm() < radius;
  if (deviation.norm() < 1.0e-5)
    deviation = Vector2D(1, 0); // does not matter, we are in the center, but avoid NaNs
  deviation.normalize();
  auto q = center + deviation * radius;
  u = q[0];
  v = q[1];
  point = Point3D(33 * q[0], 16, 33 * q[1] - 6.5);
  Vector2D perp(deviation[1], -deviation[0]);
  auto ddev = Vector3D(0, -1.2, 0) * multiplier;
  auto dperp = Vector3D(perp[0], 0, perp[1]) * multiplier;
  du = ddev * deviation[0] + dperp * perp[0];
  dv = ddev * deviation[1] + dperp * perp[1];
  return inside ? -1 : (p - q).norm();
}

double hole2(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv) {
  Point2D center(0.7, 0.5), p(u, v);
  double radius = 0.1;
  auto deviation = p - center;
  bool inside = deviation.norm() < radius;
  if (deviation.norm() < 1.0e-5)
    deviation = Vector2D(1, 0); // does not matter, we are in the center, but avoid NaNs
  deviation.normalize();
  auto q = center + deviation * radius;
  u = q[0];
  v = q[1];
  point = Point3D(33 * q[0], 16, 33 * q[1] - 6.5);
  Vector2D perp(deviation[1], -deviation[0]);
  auto ddev = Vector3D(0, -1.2, 0) * multiplier;
  auto dperp = Vector3D(perp[0], 0, perp[1]) * multiplier;
  du = ddev * deviation[0] + dperp * perp[0];
  dv = ddev * deviation[1] + dperp * perp[1];
  return inside ? -1 : (p - q).norm();
}


// Main code

int main() {
  std::vector<Boundary> boundaries = { b1, b2, b3, b4 }, holes = { hole1, hole2 };

  std::set<size_t> hole_indices;
  auto is_in_hole = [&](size_t i) { return hole_indices.find(i) != hole_indices.end(); };
  boundaries.insert(boundaries.end(), holes.begin(), holes.end());

  std::ofstream f("/tmp/sabin.obj");
  size_t index = 1;
  for (size_t i = 0; i <= resolution; ++i) {
    double u = (double)i / resolution;
    for (size_t j = 0; j <= resolution; ++j, ++index) {
      double v = (double)j / resolution;

      // Heuristic handling of holes
      bool in_hole = false;
      for (const auto &hole : holes) {
        double u1 = u, v1 = v;
        Point3D p;
        Vector3D du, dv;
        if (hole(u1, v1, p, du, dv) < 0) {
          in_hole = true;
          hole_indices.insert(index);
          Vector3D n = (du ^ dv).normalize();
          f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
          f << "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
          break;
        }
      }
      if (in_hole)
        continue;

      auto [p, n] = evalPatch(boundaries, u, v);
      f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
      f << "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
    }
  }
  for (size_t i = 0; i < resolution; ++i)
    for (size_t j = 0; j < resolution; ++j) {
      index = i * (resolution + 1) + j + 1;
      if (!is_in_hole(index) || 
          !is_in_hole(index + resolution + 2) ||
          !is_in_hole(index + 1))
        f << "f " << index << "//" << index
          << ' '  << index + resolution + 2 << "//" << index + resolution + 2
          << ' '  << index + 1 << "//" << index + 1 << std::endl;
      if (!is_in_hole(index) ||
          !is_in_hole(index + resolution + 1) ||
          !is_in_hole(index + resolution + 2))
        f << "f " << index << "//" << index
          << ' '  << index + resolution + 1 << "//" << index + resolution + 1
          << ' '  << index + resolution + 2 << "//" << index + resolution + 2 << std::endl;
    }
}
