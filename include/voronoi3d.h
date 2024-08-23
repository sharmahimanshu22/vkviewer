#ifndef VORONOI_H
#define VORONOI_H

struct Delaunay {
  std::vector<glm::dvec3> points;
  std::vector<glm::uvec4> tetrahedra; // store veretex indices // store in positive orientation
  std::vector<glm::ivec4> tetToTets;  // adjacent tetrahedra to four faces
};

Delaunay generateDelaunayTest();


#endif
