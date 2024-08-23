#include "utils.h"
#include "glm/glm.hpp"
#include <iostream>
#include <vector>

#include <algorithm>
#include <cassert>

#include <glm/gtc/type_ptr.hpp>

#include "Predicates_psm.h"
#include "voronoi3d.h"

// eight points must be passed
// maybe in some order ?

/*
  Coordinate System used
  RHS
    |z
    |
    |
    |_ _ _ _x


    minx,miny,minz  a
    maxx,miny,minz  b
    maxx,maxy,minz  c
    minx,maxy,minz  d
    minx,miny,maxz  e
    maxx,miny,maxz  f
    maxx,maxy,maxz  g
    minx,maxy,maxz  h
*/

std::vector<glm::dvec3> getBoundingBox (std::vector<glm::dvec3> points) {

  double minx, miny, minz = std::numeric_limits<double>::infinity();
  double maxx, maxy, maxz = -std::numeric_limits<double>::infinity();
  
  for (glm::dvec3 p : points) {
    if (p.x < minx) {
      minx = p.x;
    }
    if(p.y < miny) {
      miny = p.y;
    }
    if(p.z < minz) {
      minz = p.z;
    }
    if (p.x > maxx) {
      maxx = p.x;
    }
    if(p.y > maxy) {
      maxy = p.y;
    }
    if(p.z > maxz) {
      maxz = p.z;
    }
  }

  minx = minx-1;
  miny = miny-1;
  minz = minz-1;
  maxx = maxx+1;
  maxy = maxy+1;
  maxz = maxz+1;

  glm::dvec3 a(minx,miny,minz);
  glm::dvec3 b(maxx,miny,minz);
  glm::dvec3 c(maxx,maxy,minz);
  glm::dvec3 d(minx,maxy,minz);
  glm::dvec3 e(minx,miny,maxz);
  glm::dvec3 f(maxx,miny,maxz);
  glm::dvec3 g(maxx,maxy,maxz);
  glm::dvec3 h(minx,maxy,maxz);

  std::vector<glm::dvec3> bbox = {a,b,c,d,e,f,g,h};
  return bbox;
}


void constructFirstFiveTetra(Delaunay& del) {

  /*
  0124
    1327
    2476
    1745
    1724 - middle one
  */

  glm::uvec4 tet0 (0,1,3,4);
  glm::uvec4 tet1 (1,2,3,6);
  glm::uvec4 tet2 (3,4,6,7);
  glm::uvec4 tet3 (1,6,4,5);
  glm::uvec4 tet4 (1,6,3,4);


  del.tetrahedra.resize(5);
  del.tetrahedra[0] = tet0;
  del.tetrahedra[1] = tet1;
  del.tetrahedra[2] = tet2;
  del.tetrahedra[3] = tet3;
  del.tetrahedra[4] = tet4;
  
  del.tetToTets.resize(5); 
  
  del.tetToTets[0] = glm::ivec4(4,-1,-1,-1);
  del.tetToTets[1] = glm::ivec4(-1,4,-1,-1);
  del.tetToTets[2] = glm::ivec4(-1,-1,-1,4);
  del.tetToTets[3] = glm::ivec4(-1,-1,-1,4);
  del.tetToTets[4] = glm::ivec4(2,0,3,1);
    
    

}

glm::uvec3 getFacet(uint8_t i, uint32_t tetIdx, Delaunay& del) {
  //assert("i shold satisfy 1 <= 4\n", i <= 4);

  glm::uvec4 tet = del.tetrahedra[tetIdx];

  glm::uvec3 face;
  if (i == 0) {
    face.x = tet[1];
    face.y = tet[2];
    face.z = tet[3];
  }
  if (i ==1) {
    face.x = tet[0];
    face.y = tet[3];
    face.z = tet[2];
  }
  if (i == 2) {
    face.x = tet[3];
    face.y = tet[0];
    face.z = tet[1];
  }
  if (i==3) {
    face.x = tet[1];
    face.y = tet[0];
    face.z = tet[2];
  }
  return face;
}



std::vector<glm::dvec3> getPoints() {
  std::vector<glm::dvec3> points = generate3dPoints(50, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0);
  std::cout << "points size : " << points.size() << "\n"; 
  // If we want it sorted by z axis
  //std::vector<glm::vec3> sortedPoints =
  //  std::sort(points.begin(), points.end(), [](const vec_3d_t& lhs, const vec_3d_t& rhs) {
  //    return lhs.z < rhs.z;
  //  });
  return points;
}


void flip23(Delaunay& del, uint32_t tetIdx1, uint32_t tetIdx2) {

  glm::uvec4 tet1 = del.tetrahedra[tetIdx];
  glm::ivec4 adjTetIdces1 = del.tetToTets[tetIdx];

  glm::uvec4 tet2 = del.tetrahedra[tetIdx];
  glm::ivec4 adjTetIdces2 = del.tetToTets[tetIdx];

  p1, p2, a, b, c

    p1,p2,a,b
    p1,p2,b,c
    p1,p2,c,a
  
}


void flip32(Delaunay& del, uint32_t tetIdx1, uint32_t tetIdx2) {
  

}

void flip14(Delaunay& del, uint32_t tetIdx, uint32_t ptIdx) {

  glm::uvec4 tet = del.tetrahedra[tetIdx];
  glm::ivec4 adjTetIdces = del.tetToTets[tetIdx];
  
  //glm::uvec3 faceVtxIdcesOppA = getFacet(0, tetIdx, del);
  //glm::uvec3 faceVtxIdcesOppB= getFacet(1, tetIdx, del);
  //glm::uvec3 faceVtxIdcesOppC = getFacet(2, tetIdx, del);
  //glm::uvec3 faceVtxIdcesOppD = getFacet(3, tetIdx, del);

  del.tetrahedra[tetIdx] = glm::dvec4(ptIdx, tet[1], tet[2], tet[3]);
  uint32_t oppAIdx = tetIdx;
  
  del.tetrahedra.push_back(glm::dvec4(ptIdx, tet[0], tet[3], tet[2]));
  uint32_t oppBIdx = del.tetrahedra.size()-1;
  
  del.tetrahedra.push_back(glm::dvec4(ptIdx, tet[3], tet[0],tet[1]));
  uint32_t oppCIdx = del.tetrahedra.size()-1;

  del.tetrahedra.push_back(glm::dvec4(ptIdx, tet[1], tet[0], tet[2]));
  uint32_t oppDIdx = del.tetrahedra.size()-1;

  glm::vec4 l(oppAIdx, oppBIdx, oppCIdx, oppDIdx);

  del.tetToTets.resize(del.tetToTets.size() + 3);

  //del.tetToTets[oppAIdx] = glm::vec4(adjTetIdces[0], l[1], l[2], l[3]);
  //del.tetToTets[oppBIdx] = glm::vec4(l[0], adjTetIdces[1], l[2], l[3]);
  //del.tetToTets[oppCIdx] = glm::vec4(l[0], l[1],adjTetIdces[2], l[3]);
  //del.tetToTets[oppDIdx] = glm::vec4(l[0], l[1], l[2], adjTetIdces[3]);
  
  del.tetToTets[oppAIdx] = glm::vec4(adjTetIdces[0], l[1], l[2], l[3]);
  del.tetToTets[oppBIdx] = glm::vec4(adjTetIdces[1], l[0], l[3], l[2]);
  del.tetToTets[oppCIdx] = glm::vec4(adjTetIdces[2], l[3], l[0], l[1]);
  del.tetToTets[oppDIdx] = glm::vec4(adjTetIdces[3], l[1], l[0], l[2]);

  int oppATetIdx = adjTetIdces[0];
  for (int i = 0 ; i < 4 ; i++) {
    if ( del.tetToTets[oppATetIdx][i] == tetIdx) {
      del.tetToTets[oppATetIdx][i] = oppAIdx;
    }
  }
  int oppBTetIdx = adjTetIdces[1]; 
  for (int i = 0 ; i < 4 ; i++) {
    if ( del.tetToTets[oppBTetIdx][i] == tetIdx) {
      del.tetToTets[oppBTetIdx][i] = oppBIdx;
    }
  }
  int oppCTetIdx = adjTetIdces[2]; 
  for (int i = 0 ; i < 4 ; i++) {
    if ( del.tetToTets[oppCTetIdx][i] == tetIdx) {
      del.tetToTets[oppCTetIdx][i] = oppCIdx;
    }
  }
  int oppDTetIdx = adjTetIdces[3]; 
  for (int i = 0 ; i < 4 ; i++) {
    if ( del.tetToTets[oppDTetIdx][i] == tetIdx) {
      del.tetToTets[oppDTetIdx][i] = oppDIdx;
    }
  }

}

  
uint32_t locate(Delaunay del, glm::dvec3 point, int startTetIdx, bool& onFace) {

  
  if (startTetIdx < 0) {
    startTetIdx = rand() % (del.tetrahedra.size());
  }
  int tetIdx = startTetIdx;
  bool onFaceTemp;
  int res;
  uint32_t count = 0;
  while(true) {
    std::cout << tetIdx << " " << del.tetrahedra.size() << " wow\n";

    if(count > del.tetrahedra.size()) {
      throw std::runtime_error("while loop more than size of delaunay\n");
    }
    if(tetIdx < 0 || tetIdx >= del.tetrahedra.size()) {
      //std::cout << tetIdx << " : here\n";
      throw std::runtime_error("\ntetIdx is not in limits\n");
    }
    count++;
    
    glm::uvec4 tet = del.tetrahedra[tetIdx];
    double* pt = glm::value_ptr(point);
    double* vtx0= glm::value_ptr(del.points[tet[0]]);
    double* vtx1= glm::value_ptr(del.points[tet[1]]);
    double* vtx2= glm::value_ptr(del.points[tet[2]]);
    double* vtx3= glm::value_ptr(del.points[tet[3]]);
    onFaceTemp = false;

    
    res = GEO::PCK::orient_3d(pt,vtx1,vtx2,vtx3);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][0];
      continue;
    }
    if(res == 0) {
      onFaceTemp = true;
    }
    
    res = GEO::PCK::orient_3d(vtx0, pt, vtx2, vtx3);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][1];
      continue;
    }
    if(res == 0) {
      onFaceTemp = true;
    }
    
    res = GEO::PCK::orient_3d(vtx0, vtx1, pt, vtx3);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][2];
      continue;
    }
    if(res == 0) {
      onFaceTemp = true;
    }
    
    res = GEO::PCK::orient_3d(vtx0, vtx1, vtx2, pt);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][3];
      continue;
    }
    if(res == 0) {
      onFaceTemp = true;
    }

    onFace = onFaceTemp;
    break;
  }

  std::cout << tetIdx << " result of tet\n";
  return tetIdx;
  
}

std::ostream &operator<<(std::ostream &os, glm::dvec3 const &m) { 
  return os << m.x << " " << m.y << " " << m.z << "\n";
}

void insertPoint(Delaunay& del, uint32_t i) {

  bool onFace;
  uint32_t tetIdx = locate(del, del.points[i], -1, onFace);
  flip14(del, tetIdx, i);
  
}


Delaunay generateDelaunay(const std::vector<glm::dvec3>& points) {
  Delaunay del;
  del.points.clear();
  del.tetrahedra.clear();
  del.tetToTets.clear();

  std::vector<glm::dvec3> bbox = getBoundingBox(points);

  //for(auto d : bbox) {
  // std::cout << d << "\n";
  //}
  
  del.points.resize(points.size() + 8);

  std::copy(bbox.begin(), bbox.end(), del.points.begin());
  std::copy(points.begin(), points.end(), del.points.begin()+8);

  std::cout << del.points.size() <<"reached here\n";
  constructFirstFiveTetra(del);
  
  for(int i = 8; i < del.points.size(); i++) {
    insertPoint(del, i);
  }

  std::cout << del.points.size() << " " << del.tetrahedra.size() << " points and tets\n";

  return del;
}


Delaunay generateDelaunayTest() {
  std::vector<glm::dvec3> points = getPoints();
  Delaunay del = generateDelaunay(points);
  return del;
}
  




/*
std::vector<glm::dvec3>  getInfTetraCoordinates() {

  std::vector<glm::dvec3> points(4);
  points[0].x = 0.0;
  points[0].y = 0.0;
  points[0].z = std::numeric_limits<double>::infinity();
  points[1].x = 0.0;
  points[1].y = std::numeric_limits<double>::infinity();
  points[1].z = -std::numeric_limits<double>::infinity();
  points[2].x = std::numeric_limits<double>::infinity();
  points[2].y = -std::numeric_limits<double>::infinity();
  points[2].z = -std::numeric_limits<double>::infinity();
  points[3].x = -std::numeric_limits<double>::infinity();
  points[3].y = -std::numeric_limits<double>::infinity();
  points[3].z = -std::numeric_limits<double>::infinity();
  return points;

}
*/
