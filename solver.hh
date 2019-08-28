#pragma once

#include <geometry.hh>

using namespace Geometry;

using Boundary =
  double(*)(double &u, double &v, Point3D &point, Vector3D &du, Vector3D &dv);

std::pair<Point3D, Vector3D>
evalPatch(const std::vector<Boundary> &boundaries, double u, double v);
