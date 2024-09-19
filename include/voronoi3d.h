#ifndef VORONOI_H
#define VORONOI_H


#include "point3d.h"

struct Delaunay {
  std::vector<Point3d> points;
  std::vector<glm::uvec4> tetrahedra; // store veretex indices // store in positive orientation
  std::vector<glm::ivec4> tetToTets;  // adjacent tetrahedra to four faces
};

Delaunay* generateDelaunayTest();

int do_delaunay();

#endif
