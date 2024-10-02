#include "utils.h"
#include "glm/glm.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <glm/gtc/type_ptr.hpp>
#include "Predicates_psm.h"
#include "voronoi3d.h"
#include <stack>
#include "point3d.h"
#include <set>
#include "objLoader.h"
#include "vkviewer.h"

std::set<uint32_t> emptyTetras;
uint32_t lastIdx = 0;
HelloTriangleApplication app;

std::ostream &operator<<(std::ostream &os, Point3d const &m) {
  return os << m.x << " " << m.y << " " << m.z;
}

std::ostream &operator<<(std::ostream &os, glm::uvec4 const &m) {
  return os << m[0] << " " << m[1] << " " << m[2] << " " << m[3];
}

std::ostream &operator<<(std::ostream &os, glm::ivec4 const &m) {
  return os << m[0] << " " << m[1] << " " << m[2] << " " << m[3];
}


int add_point_to_render(const Point3d pt, const float size, vkview::DataForGPU& dataForGPU) {

  vkview::Vertex vertex{};

  uint32_t offset = dataForGPU.vertices.size();
  
  vertex.pos = {pt.x, pt.y, pt.z};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);
  
  vertex.pos = {pt.x-size , pt.y, pt.z};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);
  
  vertex.pos = {pt.x+size , pt.y, pt.z};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);
  
  vertex.pos = {pt.x , pt.y-size, pt.z};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);
  
  vertex.pos = {pt.x, pt.y+size, pt.z};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);
  
  vertex.pos = {pt.x, pt.y, pt.z-size};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);
  
  vertex.pos = {pt.x, pt.y, pt.z+size};
  vertex.color = {1.0,0.0,0.0};
  vertex.texCoord = {0.0f,0.0f};
  dataForGPU.vertices.push_back(vertex);

  std::vector<uint32_t> pointIndices = {offset+2,offset+4,offset+5,offset+2,offset+4,offset+6,
				offset+2,offset+3,offset+5,offset+2,offset+3,offset+6,
				offset+1,offset+4,offset+5,offset+1,offset+4,offset+6,
				offset+1,offset+3,offset+5,offset+1,offset+3,offset+6};
  
  dataForGPU.indices.insert(dataForGPU.indices.end(), pointIndices.begin(), pointIndices.end());
  return 0;
}


int add_tetras_to_render(Delaunay* del, std::vector<uint32_t> tetIdces, vkview::DataForGPU& dataForGPU) {

  
  std::unordered_map<vkview::Vertex, uint32_t> uniqueVertices{};
  //vkview::DataForGPU dataForGPU{};

  std::vector<glm::vec3> colors;
  colors.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
  colors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
  colors.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
  
  vkview::Vertex vertex{};
  Point3d pt;
  bool boundary = false;
  int j = 0;

  if (tetIdces.size() == 0) {
    for(int i = 0; i < del->tetrahedra.size(); i++) {
      tetIdces.push_back(i);
    }
  }
  
  int k = 0;
  for(int i = 0; i < tetIdces.size(); i++) {
    
    uint32_t tetIdx =  tetIdces[i];
    glm::uvec4& tet = del->tetrahedra[tetIdx];
    //for(const glm::uvec4& tet : del.tetrahedra) {
    
    uint32_t triangleIdxList[12]  = {tet[1],tet[2],tet[3],tet[0], tet[3], tet[2], tet[3],tet[0],tet[1],tet[1], tet[0], tet[2]};
    
    bool ghost = false;
    double d = 99;
    for(int i = 0; i < 12; i++) {
      pt = del->points[triangleIdxList[i]];
      if(pt.x < -d || pt.y < -d || pt.z < -d || pt.x > d || pt.y > d || pt.z > d) {
	      ghost = true;
      }
    }
    //ghost = false;
    if(ghost) {
      ghost = false;
      continue;
    }
    for(int j = 0 ; j < 12; j++) {
      pt = del->points[triangleIdxList[j]];
      vertex.pos = {pt.x , pt.y, pt.z};
      vertex.color = colors[0];
      k++;
      vertex.texCoord = {0.0f,0.0f};
      
      if (uniqueVertices.count(vertex) == 0) {
	      uniqueVertices[vertex] = static_cast<uint32_t>(dataForGPU.vertices.size());
	      dataForGPU.vertices.push_back(vertex);
      }
      dataForGPU.indices.push_back(uniqueVertices[vertex]);
    }
  }

  return 0;
}


int draw(vkview::DataForGPU& dataForGPU) {
  
  //return dataForGPU;
  try {
    app.loadScene(dataForGPU);
    app.run();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
  
}

void draw_tetras_and_point(Delaunay* del, const Point3d& p) {
  vkview::DataForGPU dataForGPU{};
  std::vector<uint32_t> emptyvector;
  add_tetras_to_render(del, emptyvector, dataForGPU);
  add_point_to_render(p, 0.05, dataForGPU);
  draw(dataForGPU);
}

void draw_tetras_and_points(Delaunay* del, std::vector<uint32_t> tetIdces, const std::vector<uint32_t>& pts) {
  vkview::DataForGPU dataForGPU{};
  add_tetras_to_render(del, tetIdces, dataForGPU);
  for(const uint32_t pix :  pts) {
    add_point_to_render(del->points[pix], 0.05, dataForGPU);
  }
  draw(dataForGPU);
}

void draw_tetras_and_point(Delaunay* del, std::vector<uint32_t> tetIdces, const Point3d& p) {
  vkview::DataForGPU dataForGPU{};
  add_tetras_to_render(del, tetIdces, dataForGPU);
  add_point_to_render(p, 0.05, dataForGPU);
  draw(dataForGPU);
}

int draw_tetra(Delaunay* del, std::vector<uint32_t> tetIdces) {
  vkview::DataForGPU dataForGPU{};
  add_tetras_to_render(del, tetIdces, dataForGPU);
  return draw(dataForGPU);
}

int draw_tetra(Delaunay* del) {
  std::vector<uint32_t> emptyvector;
  return draw_tetra(del, emptyvector);
}

std::vector<Point3d> get_points(int N, int seed) {

  /*
  std::vector<Point3d> points(4);
  
  Point3d p1 = {-1.0, -1.0, -1.0};
  Point3d p2 = {1.0, -1.0, -1.0};
  Point3d p3 = {0.0, 1.0, -1.0};
  Point3d p4 = {0.0,0.0, 1.0};

  //Point3d p5 = {0.0, -1.0, -1.0};

  points[0] = p1;
  points[1] = p2;
  points[2] = p3;
  points[3] = p4;
  //points[4] = p5;
  */
  
  std::vector<Point3d> points = generate3dPoints(N, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0, seed);

  // If we want it sorted by z axis
  //std::vector<glm::vec3> sortedPoints =
  //  std::sort(points.begin(), points.end(), [](const vec_3d_t& lhs, const vec_3d_t& rhs) {
  //    return lhs.z < rhs.z;
  //  });
  return points;
}

int orient_3d(Point3d p1, Point3d p2, Point3d p3, Point3d p4) {
  double a[3] = {p1.x, p1.y, p1.z};
  double b[3] = {p2.x, p2.y, p2.z};
  double c[3] = {p3.x, p3.y, p3.z};
  double d[3] = {p4.x, p4.y, p4.z};
  return GEO::PCK::orient_3d(a,b,c,d);
}

glm::uvec3 get_opp_face_local_idces(uint8_t i) {
  if(i==0) {
    return glm::uvec3(1,2,3);
  }
  if(i==1) {
    return glm::uvec3(0,3,2);
  }
  if(i==2) {
    return glm::uvec3(3,0,1);
  }
  if(i==3) {
    return glm::uvec3(1,0,2);
  }
  throw std::runtime_error("parameter should be between 0 and 3.\n");
}

uint32_t find_local_from_global(Delaunay* del, uint32_t tetIdx, uint32_t idx) {
  glm::uvec4 tet = del->tetrahedra[tetIdx];
  for(int i = 0; i < 4; i++) {
    if (tet[i] == idx) {
      return i;
    }
  }
  std::vector<uint32_t> ttd = {17}; 
  draw_tetras_and_point(del, ttd, del->points[idx]);
  throw std::runtime_error(" global not found in local \n");
}

uint32_t add_tetra(Delaunay* del, glm::uvec4& tet) {
  if(emptyTetras.empty()) {
    del->tetrahedra.push_back(tet);
    return del->tetrahedra.size()-1;
  } else {
    uint32_t idx = *emptyTetras.begin();
    emptyTetras.erase(emptyTetras.begin());
    del->tetrahedra[idx] = tet;
    return idx;
  }
}

// will not work for a degenerate case where two vertices have same global index
// also stupid inputs will not cause error but will lead to stupid outputs
uint32_t get_local_vtx_idx_opp_to_face(glm::uvec4 tet, uint32_t v1, uint32_t v2, uint32_t v3) {
  for(int i = 0; i < 4; i++) {
    if(tet[i] == v1 || tet[i] == v2 || tet[i] == v3) {
      continue;
    }
    return i;
  }
  throw std::runtime_error("the face supplied is not \n");
}

void printTet(glm::uvec4 tet) {
  std::cout << tet[0] << " " << tet[1] << " " << tet[2] <<  " " << tet[3] << "\n";
}

void assert_consistency(Delaunay* del) {

  assert(del->tetrahedra.size() == del->tetToTets.size());
  
  for(int tetIdx = 0; tetIdx < del->tetrahedra.size(); tetIdx++) {
    if (emptyTetras.find(tetIdx) != emptyTetras.end()) {
      continue;
    }

    glm::ivec4 adjIdces = del->tetToTets[tetIdx];
    for(int a = 0; a < 4; a++) {
      for (int b = a + 1; b < 4; b++) {
        assert (adjIdces[a] == -1 || adjIdces[a] != adjIdces[b]);
      }
    }
    for(int j = 0; j < 4; j++) {
      int adjIdx = del->tetToTets[tetIdx][j];
      assert (adjIdx != tetIdx);
      assert(emptyTetras.find(adjIdx) == emptyTetras.end());

      if(adjIdx == -1) {
	      continue;
      }
      glm::uvec4 tet = del->tetrahedra[tetIdx];
      glm::uvec4 adjTet = del->tetrahedra[adjIdx];

      std::set<uint32_t> common;
      for (int k = 0; k < 4; k++) {
        //assert(adjTet[k] != tet[j]);
	      common.insert(tet[k]);
	      common.insert(adjTet[k]);
      }
      assert(common.size() == 5);
      
    }
  }

  
}


void preprocess_points(std::vector<Point3d>& points) {
  int res;
  double pt0[3] = {points[0].x,points[0].y,points[0].z};
  double pt1[3] = {points[1].x,points[1].y,points[1].z};
  double pt2[3] = {points[102].x,points[102].y,points[102].z};
  
  res = GEO::PCK::orient_2d(pt0, pt1, pt2);
  res = 1;
  int swapIdx = 3;
  while(res == 0) {
    Point3d temp = points[swapIdx];
    points[swapIdx] = points[2];
    points[2] = temp;
    pt2[0] = points[2].x;
    pt2[1] = points[2].y;
    pt2[2] = points[2].z;
    res = GEO::PCK::orient_2d(pt0, pt1, pt2);
    temp = points[swapIdx];
    points[swapIdx] = points[2];
    points[2] = temp;

    swapIdx++;
    if(swapIdx >= points.size()) {
      throw std::runtime_error("not enough non collinear points\n");
    }
  }
  
  double pt3[3] = {points[3].x,points[3].y,points[3].z};
  
  res = GEO::PCK::orient_3d(pt0, pt1, pt2, pt3); 
  swapIdx = 4;
  while(res == 0) {
    Point3d temp = points[swapIdx];
    points[swapIdx] = points[3];
    points[3] = temp;
    pt3[0] = points[3].x;
    pt3[1] = points[3].y;
    pt3[2] = points[3].z;
    res = GEO::PCK::orient_3d(pt0, pt1, pt2, pt3);
    swapIdx++;
    if(swapIdx >= points.size()) {
      throw std::runtime_error("not enough non coplanar points\n");
    }
  }
  return;
}

void initialize_delaunay_naive(std::vector<Point3d>& points, Delaunay* del) {
  del->points.clear();
  del->tetrahedra.clear();
  del->tetToTets.clear();
  //preprocess_points(points);
  del->points.resize(points.size()+4);
  std::copy(points.begin(), points.end(), del->points.begin()+4);

  double d = 1000;
  del->points[0] = {-d, -d, -d};
  del->points[1] = {d, -d, -d};
  del->points[2] = {0.0, d, -d};
  del->points[3] = {0.0,0.0, d};

  glm::uvec4 tet0(0,1,2,3);

  del->tetrahedra.push_back(tet0);
  del->tetToTets.resize(1);
  del->tetToTets[0] = glm::ivec4(-1,-1,-1,-1);    
}

uint32_t find_local_idx_opp_to_global_face(Delaunay* del, uint32_t tetIdx, uint32_t a, uint32_t b, uint32_t c) {
  glm::uvec4 tet = del->tetrahedra[tetIdx];
  uint32_t h;
  for(int i = 0; i < 4; i++) {
    if(tet[i] != a && tet[i] != b && tet[i] != c) {
      return i;
    }
  }
  throw std::runtime_error("wrong question \n");
}

int get_orientation_wrt_tetra_with_replacement(Delaunay* del, const uint32_t tetIdx, const Point3d& p, const uint32_t replaceIdx) {

  assert(replaceIdx >= 0 && replaceIdx <=3);
  glm::uvec4 tet = del->tetrahedra[tetIdx];
  Point3d parr[4] = {del->points[tet[0]],del->points[tet[1]],del->points[tet[2]],del->points[tet[3]]};
  parr[replaceIdx] = p;
  return orient_3d(parr[0], parr[1], parr[2], parr[3]);
  
}


// more serious methods hereon
bool assert_flip32() {
  return true;
}

void flip32(Delaunay* del, uint32_t tetIdx1, uint32_t tetIdx2, uint32_t tetIdx3, uint32_t pIdxTet1, uint32_t dIdxTet2,
	    uint32_t commonVtxTet2, uint32_t kTet2, uint32_t lTet2, std::stack<uint32_t>& stack) {

  //std::cout << "performing flip 32 on " << tetIdx1 << " and " << tetIdx2 << " and " << tetIdx3 << "\n";
  



  std::vector<uint32_t> ttd = {5,0,6};
  //draw_tetra(del,ttd);
  int pIdxGlobal = del->tetrahedra[tetIdx1][pIdxTet1]; // pidxlocal should be 0
  int dIdxGlobal = del->tetrahedra[tetIdx2][dIdxTet2];
  
  glm::ivec4 tet1 = del->tetrahedra[tetIdx1];
  glm::ivec4 tet2 = del->tetrahedra[tetIdx2];
  glm::ivec4 tet3 = del->tetrahedra[tetIdx3];

  glm::ivec4 tet1AdjIdces = del->tetToTets[tetIdx1];
  glm::ivec4 tet2AdjIdces = del->tetToTets[tetIdx2];
  glm::ivec4 tet3AdjIdces = del->tetToTets[tetIdx3];

  glm::uvec4 newTet1;
  glm::uvec4 newTet2;
  glm::ivec4 newTet1AdjIdces;
  glm::ivec4 newTet2AdjIdces;
  
  uint32_t kTet1 = find_local_from_global(del, tetIdx1, del->tetrahedra[tetIdx2][kTet2]);
  uint32_t kTet3 = find_local_from_global(del, tetIdx3, del->tetrahedra[tetIdx2][kTet2]);
  uint32_t lTet1 = find_local_from_global(del, tetIdx1, del->tetrahedra[tetIdx2][lTet2]);
  uint32_t lTet3 = find_local_from_global(del, tetIdx3, del->tetrahedra[tetIdx2][lTet2]);

  int res;

  // newTet1
  res = orient_3d(del->points[pIdxGlobal],del->points[tet2[commonVtxTet2]],del->points[tet2[kTet2]],del->points[dIdxGlobal]);
  if(res >= 0) {
    newTet1 = glm::ivec4(pIdxGlobal, tet2[commonVtxTet2], tet2[kTet2], dIdxGlobal);
    newTet1AdjIdces = glm::uvec4(tet2AdjIdces[lTet2], tet3AdjIdces[lTet3], tetIdx2,  tet1AdjIdces[lTet1]);
  } else {
    newTet1 = glm::ivec4(pIdxGlobal, tet2[kTet2], tet2[commonVtxTet2], dIdxGlobal);
    newTet1AdjIdces = glm::uvec4(tet2AdjIdces[lTet2], tetIdx2, tet3AdjIdces[lTet3], tet1AdjIdces[lTet1]);
  }


  // newTet2
  res = orient_3d(del->points[pIdxGlobal],del->points[tet2[commonVtxTet2]],del->points[tet2[lTet2]],del->points[dIdxGlobal]);
  if(res >= 0) {
    newTet2 = glm::ivec4(pIdxGlobal, tet2[commonVtxTet2], tet2[lTet2], dIdxGlobal);
    newTet2AdjIdces = glm::uvec4(tet2AdjIdces[kTet2], tet3AdjIdces[kTet3], tetIdx1,  tet1AdjIdces[kTet1]);
  } else {
    newTet2 = glm::ivec4(pIdxGlobal, tet2[lTet2], tet2[commonVtxTet2], dIdxGlobal);
    newTet2AdjIdces = glm::uvec4(tet2AdjIdces[kTet2], tetIdx1, tet3AdjIdces[kTet3], tet1AdjIdces[kTet1]);
  }



  del->tetrahedra[tetIdx1] = newTet1;
  del->tetrahedra[tetIdx2] = newTet2; 
  del->tetToTets[tetIdx1] = newTet1AdjIdces;
  del->tetToTets[tetIdx2] = newTet2AdjIdces;
  

  uint32_t oppTetToKTet1 = tet1AdjIdces[kTet1];
  if (oppTetToKTet1 != -1) {
    for(int i = 0; i < 4; i++) {
      if(del->tetToTets[oppTetToKTet1][i] == tetIdx1) {
        del->tetToTets[oppTetToKTet1][i] = tetIdx2;
      } 
    }
  }

  uint32_t oppTetToLTet2 = tet2AdjIdces[lTet2];
  if(oppTetToLTet2 != -1) {
    for(int i = 0; i < 4; i++) {
      if(del->tetToTets[oppTetToLTet2][i] == tetIdx2) {
        del->tetToTets[oppTetToLTet2][i] = tetIdx1;
      }
    }
  }



  uint32_t oppTetToKTet3 = tet3AdjIdces[kTet3];
  if(oppTetToKTet3 != -1) {
    for(int i = 0; i < 4; i++) {
      if(del->tetToTets[oppTetToKTet3][i] == tetIdx3) {
        del->tetToTets[oppTetToKTet3][i] = tetIdx2;
      }
    }
  }

  uint32_t oppTetToLTet3 = tet3AdjIdces[lTet3];
  if(oppTetToLTet3 != -1) {
    for(int i = 0; i < 4; i++) {
      if(del->tetToTets[oppTetToLTet3][i] == tetIdx3) {
        del->tetToTets[oppTetToLTet3][i] = tetIdx1;
      }
    }
  }

  stack.push(tetIdx1);
  stack.push(tetIdx2);;

  emptyTetras.insert(tetIdx3);
  lastIdx = tetIdx2;



  //assert_consistency(del);
}



bool flip23_assert(Delaunay* del, uint32_t tetIdx1, uint32_t tetIdx2, uint32_t pIdxTet1, uint32_t dIdxTet2) {


  glm::uvec4 tet1 = del->tetrahedra[tetIdx1];
  glm::uvec4 tet2 = del->tetrahedra[tetIdx2];

  int res = get_orientation_wrt_tetra_with_replacement(del, tetIdx2, del->points[tet1[pIdxTet1]], dIdxTet2);
  assert(res <= 0);

  for(int i = 0; i < 4; i++) {
    if(i == dIdxTet2) {
      continue;
    }
    res = get_orientation_wrt_tetra_with_replacement(del, tetIdx2, del->points[tet1[pIdxTet1]], i);
    assert(res >= 0);
  }
  
  std::set<uint32_t> tet1commonface;
  std::set<uint32_t> tet2commonface;
  for(int i = 0; i < 4; i++) {
    if(i != pIdxTet1) {
      tet1commonface.insert(tet1[i]);
    }
    if(i != dIdxTet2) {
      tet2commonface.insert(tet2[i]);
    }
  }

  assert(tet1commonface == tet2commonface);

  return true;
}

// tet1 should have p
// tet2 should have d
// localIdx1 should always be 0 for now
void flip23(Delaunay* del, uint32_t tetIdx1, uint32_t tetIdx2, uint32_t pIdxTet1, uint32_t dIdxTet2, std::stack<uint32_t>& stack) {

  //std::cout << "performing flip 23 on " << tetIdx1 << " and " << tetIdx2 << " \n";
  assert(pIdxTet1 == 0);

  glm::uvec4 tet1 = del->tetrahedra[tetIdx1];
  glm::uvec4 tet2 = del->tetrahedra[tetIdx2];


  
  // check if the two tets provided are indeed applicable
  
  
  flip23_assert(del, tetIdx1, tetIdx2, pIdxTet1, dIdxTet2);
  
  uint32_t idx1_in_tet2 = find_local_from_global(del, tetIdx2, tet1[1]);
  uint32_t idx2_in_tet2 = find_local_from_global(del, tetIdx2, tet1[2]);
  uint32_t idx3_in_tet2 = find_local_from_global(del, tetIdx2, tet1[3]);

  glm::ivec4 tet1Adj = del->tetToTets[tetIdx1];
  glm::ivec4 tet2Adj = del->tetToTets[tetIdx2];
  
  int pIdx = tet1[pIdxTet1];
  int dIdx = tet2[dIdxTet2];


  // we take it on trust that the method is called for right configuration.
  glm::uvec4 newTet1 (pIdx, tet1[1], tet1[2], dIdx);
  glm::uvec4 newTet2 (pIdx, tet1[2], tet1[3], dIdx);
  glm::uvec4 newTet3 (pIdx, tet1[3], tet1[1], dIdx);

  del->tetrahedra[tetIdx1] = newTet1;
  del->tetrahedra[tetIdx2] = newTet2;
  uint32_t idx1 = tetIdx1;
  uint32_t idx2 = tetIdx2;
  uint32_t idx3 = add_tetra(del, newTet3);

  
  // updating adjacency info of new tetrahedra
  del->tetToTets[idx1] = glm::ivec4(tet2Adj[idx3_in_tet2], idx2, idx3, tet1Adj[3] );
  del->tetToTets[idx2] = glm::ivec4(tet2Adj[idx1_in_tet2], idx3, idx1, tet1Adj[1] );
  glm::ivec4 adjNewTet3(tet2Adj[idx2_in_tet2], idx1, idx2, tet1Adj[2] );
  if (idx3 == del->tetToTets.size()) {
    del->tetToTets.resize(del->tetToTets.size() + 1);
  }
  del->tetToTets[idx3] = adjNewTet3;

 
  // updating adjency info for tetras adjacent to first tetra  
  uint32_t tet1Adj_1_idx = tet1Adj[1];
  if(tet1Adj_1_idx != -1) {
    for(int j = 0; j < 4; j++) {
      if(del->tetToTets[tet1Adj_1_idx][j] == tetIdx1) {
	      del->tetToTets[tet1Adj_1_idx][j] = idx2;
      }
    }
  }
  
  uint32_t tet1Adj_2_idx = tet1Adj[2];
  if(tet1Adj_2_idx != -1) { 
    for(int j = 0; j < 4; j++) {
      if(del->tetToTets[tet1Adj_2_idx][j] == tetIdx1) {
	      del->tetToTets[tet1Adj_2_idx][j] = idx3;
      }
    }
  }
  

  uint32_t tet1Adj_3_idx = tet1Adj[3];
  if(tet1Adj_3_idx != -1) { 
    for(int j = 0; j < 4; j++) {
      if(del->tetToTets[tet1Adj_3_idx][j] == tetIdx1) {
	      del->tetToTets[tet1Adj_3_idx][j] = idx1;
      }
    }
  }
  
  // updating adjacency info for tetras adjacent to second tetra
  uint32_t tet2Adj_1_idx = tet2Adj[idx1_in_tet2];
  if(tet2Adj_1_idx != -1) { 
    for(int j = 0; j < 4; j++) {
      if(del->tetToTets[tet2Adj_1_idx][j] == tetIdx2) {
	      del->tetToTets[tet2Adj_1_idx][j] = idx2;
      }
    }
  }

  uint32_t tet2Adj_2_idx = tet2Adj[idx2_in_tet2];
  if(tet2Adj_2_idx != -1) { 
    for(int j = 0; j < 4; j++) {
      if(del->tetToTets[tet2Adj_2_idx][j] == tetIdx2) {
	      del->tetToTets[tet2Adj_2_idx][j] = idx3;
      }
    }
  }

  uint32_t tet2Adj_3_idx = tet2Adj[idx3_in_tet2];
  if(tet2Adj_3_idx != -1) { 
    for(int j = 0; j < 4; j++) {
      if(del->tetToTets[tet2Adj_3_idx][j] == tetIdx2) {
	      del->tetToTets[tet2Adj_3_idx][j] = idx1;
      }
    }
  }
  
  stack.push(idx1);
  stack.push(idx2);
  stack.push(idx3);


  //assert_consistency(del);
}





void flip44(Delaunay* del, uint32_t tetIdx1, uint32_t tetIdx2, uint32_t tetIdx3, uint32_t tetIdx4, uint32_t pIdxTet1,
	    uint32_t dIdxTet2, uint32_t commonVtxTet2, uint32_t kTet2, uint32_t lTet2, std::stack<uint32_t>& stack) {


  //std::cout << "performing flip 44 on " << tetIdx1 << " and " << tetIdx2 << " and " << tetIdx3 << " and " << tetIdx4 << "\n";

  // kl is the common edge
  glm::uvec4 tet1 = del->tetrahedra[tetIdx1];
  glm::uvec4 tet2 = del->tetrahedra[tetIdx2];
  glm::uvec4 tet3 = del->tetrahedra[tetIdx3];
  glm::uvec4 tet4 = del->tetrahedra[tetIdx4];

  glm::ivec4 oldTet1Adj = del->tetToTets[tetIdx1];
  glm::ivec4 oldTet2Adj = del->tetToTets[tetIdx2];
  glm::ivec4 oldTet3Adj = del->tetToTets[tetIdx3];
  glm::ivec4 oldTet4Adj = del->tetToTets[tetIdx4];

  uint32_t pIdxGlobal = tet1[pIdxTet1];
  uint32_t dIdxGlobal = tet2[dIdxTet2];
  uint32_t kIdxGlobal = tet2[kTet2];
  uint32_t lIdxGlobal = tet2[lTet2];
  uint32_t commonVtxTet2Global = tet2[commonVtxTet2];
  uint32_t commonVtxTet4 = find_local_idx_opp_to_global_face(del, tetIdx4, kIdxGlobal, lIdxGlobal, dIdxGlobal);
  uint32_t commonVtxTet4Global = tet4[commonVtxTet4];

  uint32_t kTet1 = find_local_from_global(del, tetIdx1, kIdxGlobal);
  uint32_t lTet1 = find_local_from_global(del, tetIdx1, lIdxGlobal);
  uint32_t kTet3 = find_local_from_global(del, tetIdx3, kIdxGlobal);
  uint32_t lTet3 = find_local_from_global(del, tetIdx3, lIdxGlobal);
  uint32_t kTet4 = find_local_from_global(del, tetIdx4, kIdxGlobal);
  uint32_t lTet4 = find_local_from_global(del, tetIdx4, lIdxGlobal);

  glm::uvec4 newTet1;
  glm::uvec4 newTet2;
  glm::uvec4 newTet3;
  glm::uvec4 newTet4;

  glm::ivec4 newTet1Adjacent;
  glm::ivec4 newTet2Adjacent;
  glm::ivec4 newTet3Adjacent;
  glm::ivec4 newTet4Adjacent;
    
  int res;  
  res = orient_3d(del->points[pIdxGlobal],del->points[commonVtxTet2Global],del->points[kIdxGlobal],del->points[dIdxGlobal]);
  if(res >= 0) {
    newTet1 = glm::uvec4(pIdxGlobal, commonVtxTet2Global, kIdxGlobal, dIdxGlobal);
    newTet1Adjacent = glm::ivec4(oldTet2Adj[lTet2], tetIdx3, tetIdx2, oldTet1Adj[lTet1]);
  } else {
    newTet1 = glm::uvec4(pIdxGlobal,  kIdxGlobal, commonVtxTet2Global, dIdxGlobal);
    newTet1Adjacent = glm::ivec4(oldTet2Adj[lTet2], tetIdx2, tetIdx3, oldTet1Adj[lTet1]);
  }

  res = orient_3d(del->points[pIdxGlobal],del->points[commonVtxTet2Global],del->points[lIdxGlobal],del->points[dIdxGlobal]);
  if(res >= 0) {
    newTet2 = glm::uvec4(pIdxGlobal, commonVtxTet2Global, lIdxGlobal, dIdxGlobal);
    newTet2Adjacent = glm::ivec4(oldTet2Adj[kTet2], tetIdx4, tetIdx1, oldTet1Adj[kTet1]);
  } else {
    newTet2 = glm::uvec4(pIdxGlobal, lIdxGlobal, commonVtxTet2Global, dIdxGlobal);
    newTet2Adjacent = glm::ivec4(oldTet2Adj[kTet2], tetIdx1, tetIdx4, oldTet1Adj[kTet1]);
  }

  res = orient_3d(del->points[pIdxGlobal],del->points[commonVtxTet4Global],del->points[kIdxGlobal],del->points[dIdxGlobal]);
  if(res >= 0) {
    newTet3 = glm::uvec4(pIdxGlobal, commonVtxTet4Global, kIdxGlobal, dIdxGlobal);
    newTet3Adjacent = glm::ivec4(oldTet4Adj[lTet4], tetIdx1, tetIdx4, oldTet3Adj[lTet3]);
  } else {
    newTet3 = glm::uvec4(pIdxGlobal, kIdxGlobal, commonVtxTet4Global, dIdxGlobal);
    newTet3Adjacent = glm::ivec4(oldTet4Adj[lTet4], tetIdx4, tetIdx1, oldTet3Adj[lTet3]);
  }

  res = orient_3d(del->points[pIdxGlobal],del->points[commonVtxTet4Global],del->points[lIdxGlobal],del->points[dIdxGlobal]);
  if(res >= 0) {
    newTet4 = glm::uvec4(pIdxGlobal, commonVtxTet4Global, lIdxGlobal, dIdxGlobal);
    newTet4Adjacent = glm::ivec4(oldTet4Adj[kTet4], tetIdx2, tetIdx3, oldTet3Adj[kTet3]);
  } else {
    newTet4 = glm::uvec4(pIdxGlobal, lIdxGlobal, commonVtxTet4Global, dIdxGlobal);
    newTet4Adjacent = glm::ivec4(oldTet4Adj[kTet4], tetIdx3, tetIdx2, oldTet3Adj[kTet3]);
  }

  del->tetrahedra[tetIdx1] = newTet1;
  del->tetrahedra[tetIdx2] = newTet2;
  del->tetrahedra[tetIdx3] = newTet3;
  del->tetrahedra[tetIdx4] = newTet4;
  del->tetToTets[tetIdx1] = newTet1Adjacent;
  del->tetToTets[tetIdx2] = newTet2Adjacent;
  del->tetToTets[tetIdx3] = newTet3Adjacent;
  del->tetToTets[tetIdx4] = newTet4Adjacent;

  int tetIdxOppToTet1K = oldTet1Adj[kTet1];
  if(tetIdxOppToTet1K != -1) {
    glm::ivec4 tetIdxOppToTet1KAdj = del->tetToTets[tetIdxOppToTet1K];
    for(int i = 0; i < 4; i++) {
      if (tetIdxOppToTet1KAdj[i] == tetIdx1) {
        del->tetToTets[tetIdxOppToTet1K][i] = tetIdx2;
        break;
      }
    }
  }
  
  int tetIdxOppToTet2L = oldTet2Adj[lTet2];
  if(tetIdxOppToTet2L != -1) {
    glm::ivec4 tetIdxOppToTet2LAdj = del->tetToTets[tetIdxOppToTet2L];
    for(int i = 0; i < 4; i++) {
      if (tetIdxOppToTet2LAdj[i] == tetIdx2) {
        del->tetToTets[tetIdxOppToTet2L][i] = tetIdx1;
        break;
      }
    }
  }
  
  int tetIdxOppToTet3K = oldTet3Adj[kTet3];
  if(tetIdxOppToTet3K != -1) {
    glm::ivec4 tetIdxOppToTet3KAdj = del->tetToTets[tetIdxOppToTet3K];
    for(int i = 0; i < 4; i++) {
      if (tetIdxOppToTet3KAdj[i] == tetIdx3) {
        del->tetToTets[tetIdxOppToTet3K][i] = tetIdx4;
        break;
      }
    }
  }
  
  int tetIdxOppToTet4L = oldTet4Adj[lTet4];
  if(tetIdxOppToTet4L != -1) {
    glm::ivec4 tetIdxOppToTet4LAdj = del->tetToTets[tetIdxOppToTet4L];
    for(int i = 0; i < 4; i++) {
      if (tetIdxOppToTet4LAdj[i] == tetIdx4) {
        del->tetToTets[tetIdxOppToTet4L][i] = tetIdx3;
        break;
      }
    }
  }

  stack.push(tetIdx1);
  stack.push(tetIdx2);
  stack.push(tetIdx3);
  stack.push(tetIdx4);

  //assert_consistency(del);

}

bool flip44_check_and_call(Delaunay* del, uint32_t tet1Idx, uint32_t tet2Idx, uint32_t pIdxTet1, uint32_t dIdxTet2, uint32_t commonVtxIdxTet2, 
    uint32_t kTet2, uint32_t lTet2, std::stack<uint32_t>& stack) {
  
  glm::uvec4 tet1 = del->tetrahedra[tet1Idx];
  glm::uvec4 tet2 = del->tetrahedra[tet2Idx];
  uint32_t commonVtxIdxTet1 = find_local_from_global(del, tet1Idx, tet2[commonVtxIdxTet2]);
  uint32_t tet4Idx = del->tetToTets[tet2Idx][commonVtxIdxTet2];
  uint32_t tet3Idx = del->tetToTets[tet1Idx][commonVtxIdxTet1];

  uint32_t pIdxTet3 = find_local_from_global(del, tet3Idx, tet1[pIdxTet1]);

  uint32_t adjToTet3AtP = del->tetToTets[tet3Idx][pIdxTet3];

  if (adjToTet3AtP == tet4Idx) {
    flip44(del, tet1Idx, tet2Idx, tet3Idx, tet4Idx, pIdxTet1, dIdxTet2, commonVtxIdxTet2, kTet2, lTet2, stack);
    return true;
  }
  else {
    return false;
  }
}

void find_case_and_resolve(Delaunay* del, const uint32_t tet1Idx, const uint32_t tet2Idx, uint32_t idx1, uint32_t idx2, std::stack<uint32_t>& stack) {
  
  assert(idx1 == 0);
  glm::uvec4 tet1 = del->tetrahedra[tet1Idx];
  Point3d p = del->points[tet1[idx1]];
  glm::uvec4 tet2 = del->tetrahedra[tet2Idx];
  int ortn[4];

  ortn[0] = orient_3d(p, del->points[tet2[1]], del->points[tet2[2]], del->points[tet2[3]]);
  ortn[1] = orient_3d(del->points[tet2[0]], p, del->points[tet2[2]], del->points[tet2[3]]);
  ortn[2] = orient_3d(del->points[tet2[0]], del->points[tet2[1]], p, del->points[tet2[3]]);
  ortn[3] = orient_3d(del->points[tet2[0]], del->points[tet2[1]], del->points[tet2[2]], p);
  
  assert(ortn[idx2] == 0 || ortn[idx2] == -1);

  glm::uvec3 jkl = get_opp_face_local_idces(idx2);
  uint32_t i = idx2;
  uint32_t j = jkl[0];
  uint32_t k = jkl[1];
  uint32_t l = jkl[2];

  uint32_t jTet1 = find_local_from_global(del, tet1Idx, tet2[j]);
  uint32_t kTet1 = find_local_from_global(del, tet1Idx, tet2[k]);
  uint32_t lTet1 = find_local_from_global(del, tet1Idx, tet2[l]);


  //ortn[i]  <= 0 and all other orientation is positive. one face visible
  // case 1
  if (ortn[i] <= 0 && ortn[j] == 1 && ortn[k] == 1 && ortn[l] == 1) {
    // only one face is visible that is ith face
    flip23(del, tet1Idx, tet2Idx, idx1, idx2, stack);
    return;
  }

  

  //ortn[i]  <= 0 and only one other orientation is negative. 2 faces visible
  // case 2
  if(ortn[i] <= 0 && ortn[j] == -1 && ortn[k] == 1 && ortn[l] == 1) {
    uint32_t tet2Adj_j = del->tetToTets[tet2Idx][j];
    uint32_t tet1Adj_j = del->tetToTets[tet1Idx][jTet1];
    if(tet1Adj_j != -1 && tet2Adj_j == tet1Adj_j ) {
      flip32(del, tet1Idx, tet2Idx, tet1Adj_j, idx1, idx2, j, k, l, stack);
      return;
    }
    //i and j are visible
  }
  if(ortn[i] <= 0 && ortn[j] == 1 && ortn[k] == -1 && ortn[l] == 1) {
    uint32_t tet2Adj_k = del->tetToTets[tet2Idx][k];
    uint32_t tet1Adj_k = del->tetToTets[tet1Idx][kTet1];
    if(tet1Adj_k != -1 && tet2Adj_k == tet1Adj_k ) {
      flip32(del, tet1Idx, tet2Idx, tet1Adj_k, idx1, idx2, k, j, l, stack);
      return;
    }
  }
  if(ortn[i] <= 0 && ortn[j] == 1 && ortn[k] == 1 && ortn[l] == -1) {
    uint32_t tet2Adj_l = del->tetToTets[tet2Idx][l];
    uint32_t tet1Adj_l = del->tetToTets[tet1Idx][lTet1];
    if(tet1Adj_l != -1 && tet2Adj_l == tet1Adj_l ) {
      flip32(del, tet1Idx, tet2Idx, tet1Adj_l, idx1, idx2, l, j, k, stack);
      return;
    }
  }
  /*
  // three faces visible // the paper mentions not to do anything in this case;
  if(ortn[i] == -1 && ortn[j] == -1 && ortn[k] == -1 && ortn[l] == 1) {
    //i and j are visible
  }
  if(ortn[i] == -1 && ortn[j] == 1 && ortn[k] == -1 && ortn[l] == -1) {
    //i and k are visible
  }
  if(ortn[i] == -1 && ortn[j] == -1 && ortn[k] == 1 && ortn[l] == -1) {
    //i and l are visible
  }
  */

  // degenerate cases

  
  if(ortn[i] == -1 && (ortn[j] == 0 || ortn[k] == 0 || ortn[l] == 0)) {

    // case 3
    if(ortn[j] == 0 && (ortn[k] != 0 && ortn[l] != 0) ) {
      flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, j, k, l, stack);
      return;
    }
    // case 3
    if(ortn[k] == 0 && (ortn[j] != 0 && ortn[l] != 0) ) {
      flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, k, j, l, stack);
      return;
    }
    //case 3
    if(ortn[l] == 0 && (ortn[j] != 0 && ortn[k] != 0) ) {
      flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, l, j, k, stack);
      return;
    }

    // umm, what to do here ?
    if(ortn[j] == 0 && ortn[k] == 0 && ortn[l] != 0 ) {
      bool resflip44;
      resflip44 = flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, j, k, l, stack);
      if(!resflip44) {
        flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, k, j, l, stack);
      }
      return;
    }
    if(ortn[k] == 0 && ortn[l] == 0 && ortn[j] != 0)  {
      bool resflip44;
      resflip44 = flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, k, j, l, stack);
      if(!resflip44) {
        flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, l, j, k, stack);
      }
      return;
    }
    if(ortn[j] == 0 && ortn[l] == 0 && ortn[k] != 0)  {
      bool resflip44;
      resflip44 = flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, j, k, l, stack);
      if(!resflip44) {
        flip44_check_and_call(del, tet1Idx, tet2Idx, idx1, idx2, l, j, k, stack);
      }
      return;
    }

    // this is a funny case. not possible
    if(ortn[j] == 0 && ortn[k] == 0 && ortn[l] == 0)  {
      // duplicate vertex. should have been deleted in pre-processing.
    }

  }

  if(ortn[i] == 0 && (ortn[j] == 0 || ortn[k] == 0 || ortn[l] == 0)) {
    // case 4
    if(ortn[j] == 0 && (ortn[k] != 0 && ortn[l] != 0) ) {
      uint32_t tet2Adj_j = del->tetToTets[tet2Idx][j];
      uint32_t tet1Adj_j = del->tetToTets[tet1Idx][jTet1];
      if(tet1Adj_j != -1 && tet2Adj_j == tet1Adj_j ) {
        flip32(del, tet1Idx, tet2Idx, tet1Adj_j, idx1, idx2, j, k, l, stack);
      } else {
        flip23(del, tet1Idx, tet2Idx, idx1, idx2, stack);
      }
      return;
    }
    //case 4
    if(ortn[k] == 0 && (ortn[j] != 0 && ortn[l] != 0) ) {
      uint32_t tet2Adj_k = del->tetToTets[tet2Idx][k];
      uint32_t tet1Adj_k = del->tetToTets[tet1Idx][kTet1];
      if(tet1Adj_k != -1 && tet2Adj_k == tet1Adj_k ) {
        flip32(del, tet1Idx, tet2Idx, tet1Adj_k, idx1, idx2, k, j, l, stack);
      } else {
        flip23(del, tet1Idx, tet2Idx, idx1, idx2, stack);
      }
      return;
    }
    //case 4
    if(ortn[l] == 0 && (ortn[j] != 0 && ortn[k] != 0) ) {
      uint32_t tet2Adj_l = del->tetToTets[tet2Idx][l];
      uint32_t tet1Adj_l = del->tetToTets[tet1Idx][lTet1];
      if(tet1Adj_l != -1 && tet2Adj_l == tet1Adj_l ) {
        flip32(del, tet1Idx, tet2Idx, tet1Adj_l, idx1, idx2, l, j, k, stack);
      } else 
        flip23(del, tet1Idx, tet2Idx, idx1, idx2, stack);
      }
      return;
    }

    // new point is common with one of vtx
    if(ortn[j] == 0 && ortn[k] == 0 && ortn[l] != 0 ) {
      // duplicate vertex. should have been deleted in pre-processing.
    }
    if(ortn[k] == 0 && ortn[l] == 0 && ortn[j] != 0)  {
      // duplicate vertex. should have been deleted in pre-processing.
    }
    if(ortn[j] == 0 && ortn[l] == 0 && ortn[k] != 0)  {
      // duplicate vertex. should have been deleted in pre-processing.
    }

    // this is a funny case. not possible
    if(ortn[j] == 0 && ortn[k] == 0 && ortn[l] == 0)  {
      // duplicate vertex. should have been deleted in pre-processing. this is actually for a point sized tetrahedra.
    }
  }


// return indices of four tetrahedra added
void flip14(Delaunay* del, uint32_t tetIdx, uint32_t ptIdx, std::stack<uint32_t>& stack, uint32_t& degenerateCase) {

  //std::cout << "performing flip 14 at " << tetIdx << "  Mesh size is : " << del->tetrahedra.size() << "\n";
  
  glm::uvec4 tet = del->tetrahedra[tetIdx];
  glm::ivec4 adjTetIdces = del->tetToTets[tetIdx];
  
  del->tetrahedra[tetIdx] = glm::uvec4(ptIdx, tet[1], tet[2], tet[3]);
  uint32_t oppAIdx = tetIdx;
  stack.push(oppAIdx);

  del->tetrahedra.push_back(glm::uvec4(ptIdx, tet[0], tet[3], tet[2]));
  uint32_t oppBIdx = del->tetrahedra.size()-1;
  stack.push(oppBIdx);
			    
  del->tetrahedra.push_back(glm::uvec4(ptIdx, tet[3], tet[0], tet[1]));
  uint32_t oppCIdx = del->tetrahedra.size()-1;
  stack.push(oppCIdx);

  del->tetrahedra.push_back(glm::uvec4(ptIdx, tet[1], tet[0], tet[2]));
  uint32_t oppDIdx = del->tetrahedra.size()-1;
  stack.push(oppDIdx);
 
  glm::ivec4 l(oppAIdx, oppBIdx, oppCIdx, oppDIdx);

  del->tetToTets.resize(del->tetToTets.size() + 3);
  
  del->tetToTets[oppAIdx] = glm::ivec4(adjTetIdces[0], l[1], l[2], l[3]);
  del->tetToTets[oppBIdx] = glm::ivec4(adjTetIdces[1], l[0], l[3], l[2]);
  del->tetToTets[oppCIdx] = glm::ivec4(adjTetIdces[2], l[3], l[0], l[1]);
  del->tetToTets[oppDIdx] = glm::ivec4(adjTetIdces[3], l[1], l[0], l[2]);


  
  int oppATetIdx = adjTetIdces[0];
  if(oppATetIdx != -1) {
    for (int i = 0 ; i < 4 ; i++) {
      if ( del->tetToTets[oppATetIdx][i] == tetIdx) {
	del->tetToTets[oppATetIdx][i] = oppAIdx;
      }
    }
  }

  int oppBTetIdx = adjTetIdces[1];
  if (oppBTetIdx != -1) {
    for (int i = 0 ; i < 4 ; i++) {
      if ( del->tetToTets[oppBTetIdx][i] == tetIdx) {
	del->tetToTets[oppBTetIdx][i] = oppBIdx;
      }
    }
  }
  
  int oppCTetIdx = adjTetIdces[2];
  if(oppCTetIdx != -1) {
    for (int i = 0 ; i < 4 ; i++) {
      if ( del->tetToTets[oppCTetIdx][i] == tetIdx) {
	del->tetToTets[oppCTetIdx][i] = oppCIdx;
      }
    }
  }
  
  
  int oppDTetIdx = adjTetIdces[3];
  if(oppDTetIdx != -1) {
    for (int i = 0 ; i < 4 ; i++) {
      if ( del->tetToTets[oppDTetIdx][i] == tetIdx) {
	del->tetToTets[oppDTetIdx][i] = oppDIdx;
      }
    }
  }
  

  lastIdx = oppAIdx;

   //std::cout << "performed flip 14 " << del->tetrahedra.size() << "\n";

  //assert_consistency(del);
  
}


int in_sphere_3d_SOS(Point3d& p0, Point3d& p1, Point3d& p2, Point3d& p3, Point3d& z) {
  
  double p0arr[3] = {p0.x, p0.y, p0.z}; // this is the newly added point p
  double p1arr[3] = {p1.x, p1.y, p1.z};
  double p2arr[3] = {p2.x, p2.y, p2.z};
  double p3arr[3] = {p3.x, p3.y, p3.z};
  double zarr[3] = {z.x, z.y, z.z};
  
  int res = GEO::PCK::in_sphere_3d_SOS(p0arr, p1arr, p2arr, p3arr, zarr);
  return res;
}

int in_sphere_3d_Inexact(Point3d& p0, Point3d& p1, Point3d& p2, Point3d& p3, Point3d& z) {
  
  double p0arr[3] = {p0.x, p0.y, p0.z}; // this is the newly added point p
  double p1arr[3] = {p1.x, p1.y, p1.z};
  double p2arr[3] = {p2.x, p2.y, p2.z};
  double p3arr[3] = {p3.x, p3.y, p3.z};
  double zarr[3] = {z.x, z.y, z.z};
  
  int res = GEO::PCK::in_sphere_3d_Inexact(p0arr, p1arr, p2arr, p3arr, zarr);
  return res;
}


void flip_if_applicable_wrt_first_vertex(Delaunay* del, int tetIdx, std::stack<uint32_t>& stack, uint32_t& degenerateCase) {

  glm::uvec4 tet = del->tetrahedra[tetIdx];
  uint32_t adjTetIdx = del->tetToTets[tetIdx][0];

  
  if(adjTetIdx == -1) {
    return;
  }
  
  glm::uvec4 adjTet = del->tetrahedra[adjTetIdx];

  uint32_t adjTetOppVtxLocalIdx = get_local_vtx_idx_opp_to_face(adjTet, tet[1], tet[2], tet[3]);
  assert(adjTetOppVtxLocalIdx >= 0 && adjTetOppVtxLocalIdx <= 3);

  
  Point3d vtx0 = del->points[tet[0]];
  Point3d vtx1 = del->points[tet[1]];
  Point3d vtx2 = del->points[tet[2]];
  Point3d vtx3 = del->points[tet[3]];
  Point3d z = del->points[adjTet[adjTetOppVtxLocalIdx]];

  int res = in_sphere_3d_SOS(vtx0, vtx1, vtx2, vtx3, z);

  if (res > 0 ) {
    find_case_and_resolve(del, tetIdx, adjTetIdx, 0, adjTetOppVtxLocalIdx, stack);
  }
}

void assert_delaunay(Delaunay* del, uint32_t ptSize) {

  for(uint32_t tetIdx = 0; tetIdx < del->tetrahedra.size(); tetIdx++) {
    if(emptyTetras.find(tetIdx) != emptyTetras.end()) {
      continue;
    }
    glm::uvec4 tet = del->tetrahedra[tetIdx];
    int res = orient_3d(del->points[tet[0]], del->points[tet[1]], del->points[tet[2]], del->points[tet[3]]);
    assert (res == 1);  
  }

  for(uint32_t tetIdx = 0; tetIdx < del->tetrahedra.size(); tetIdx++) {
    if(emptyTetras.find(tetIdx) != emptyTetras.end()) {
      continue;
    }
    glm::uvec4 tet = del->tetrahedra[tetIdx];
    for(uint32_t pIdx = 0; pIdx <= ptSize; pIdx ++) {
      if (pIdx == tet[0] || pIdx == tet[1] || pIdx == tet[2] || pIdx == tet[3]) {
        continue;
      }
      int res = in_sphere_3d_Inexact(del->points[tet[0]], del->points[tet[1]], del->points[tet[2]], del->points[tet[3]], del->points[pIdx]);

      assert(res != 1);
    }
  }
}

uint32_t locate(Delaunay* delptr, const uint32_t pointIdx, int startTetIdx, uint32_t& degenerateCase) {
  
  Delaunay del = *delptr;

  Point3d point = del.points[pointIdx];

  // make start idex to be NOT a ghost tetra
  if (startTetIdx < 0) {
    startTetIdx = lastIdx;
    //startTetIdx = rand() % (del.tetrahedra.size());
  }
  
  uint32_t tetIdx = startTetIdx;
  int res;
  uint32_t count = 0;

  std::vector<uint32_t> tetIdxStack;

  while(true) {
    tetIdxStack.push_back(tetIdx);
    if(count >= del.tetrahedra.size()) {
      throw std::runtime_error("while loop more than size of delaunay\n");
    }
    if(tetIdx < 0 || tetIdx >= del.tetrahedra.size()) {
      throw std::runtime_error("\ntetIdx is not in limits\n");
    }
    count++;
    
    glm::ivec4 tet = del.tetrahedra[tetIdx];
    //assert_not_ghost(tet);
    double pt[3] = {point.x,point.y,point.z};
    double vtx0[3] = {del.points[tet[0]].x, del.points[tet[0]].y, del.points[tet[0]].z};
    double vtx1[3] = {del.points[tet[1]].x, del.points[tet[1]].y, del.points[tet[1]].z};
    double vtx2[3] = {del.points[tet[2]].x, del.points[tet[2]].y, del.points[tet[2]].z};
    double vtx3[3] = {del.points[tet[3]].x, del.points[tet[3]].y, del.points[tet[3]].z};
    degenerateCase = 0;
    res = GEO::PCK::orient_3d(pt,vtx1,vtx2,vtx3);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][0];
      continue;
    }
    if(res == 0) {
      degenerateCase += 1;
    }
    
    res = GEO::PCK::orient_3d(vtx0, pt, vtx2, vtx3);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][1];
      continue;
    }
    if(res == 0) {
      degenerateCase += 1;
    }
    
    res = GEO::PCK::orient_3d(vtx0, vtx1, pt, vtx3);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][2];
      continue;
    }
    if(res == 0) {
      degenerateCase += 1;
    }
    
    res = GEO::PCK::orient_3d(vtx0, vtx1, vtx2, pt);
    if(res < 0) {
      tetIdx = del.tetToTets[tetIdx][3];
      continue;
    }
    if(res == 0) {
      degenerateCase += 1;
    }
    break;
  }

  return tetIdx;
}


bool assert_inside_tetra(Delaunay* del, uint32_t tetIdx, Point3d& p) {
  glm::uvec4 tet = del->tetrahedra[tetIdx];
  Point3d p0 = del->points[tet[0]];
  Point3d p1 = del->points[tet[1]];
  Point3d p2 = del->points[tet[2]];
  Point3d p3 = del->points[tet[3]];

  assert (orient_3d(p0,p1,p2,p3) >= 0);
  assert (orient_3d(p,p1,p2,p3) >= 0);
  assert (orient_3d(p0,p,p2,p3) >= 0);
  assert (orient_3d(p0,p1,p,p3) >= 0);
  assert (orient_3d(p0,p1,p2,p) >= 0);

  return true;
}

void insert_point(Delaunay* del, uint32_t ptIdx, std::stack<uint32_t>& stack) {
  uint32_t degenerateCase;

  uint32_t tetIdx = locate(del, ptIdx, -1, degenerateCase);
  
  assert_inside_tetra(del, tetIdx, del->points[ptIdx]);
  
  assert(tetIdx >= 0 && tetIdx < del->tetrahedra.size());

  flip14(del, tetIdx, ptIdx, stack, degenerateCase);
  while(!stack.empty()) {
    int tetIdx = stack.top();
    stack.pop();
    if (emptyTetras.find(tetIdx) != emptyTetras.end()) {
      continue;
    }
    flip_if_applicable_wrt_first_vertex(del, tetIdx, stack, degenerateCase);
  }
}

Delaunay* generateDelaunay(std::vector<Point3d>  points) {
  Delaunay* del = new Delaunay();
  initialize_delaunay_naive(points, del);
  std::stack<uint32_t> stack;
  emptyTetras.clear();
  lastIdx = 0;
  for(int i = 4; i < del->points.size(); i++) { 
    if (i%1000 == 0) {
      std::cout << "inserting : " << i << "\n";
    }
    insert_point(del, i, stack);
    //assert_delaunay(del, i);
  }
  //assert_delaunay(del, del->points.size()-1);
  return del;
}

Delaunay* generateDelaunayTest() {
  Point3d p1 = {0.0,-1.0, 0.0};
  Point3d p2 = {0.0,1.0, 0.0};
  Point3d p3 = {-1.0,0.0, 0.0};
  Point3d p4 = {1.0,0.0, 0.0};
  Point3d p5 = {0.0,0.0, -1.0};
  Point3d p6 = {0.0,0.0, 1.0};
  //std::vector<Point3d> points = {p1,p2,p3,p4,p5,p6};
  //std::vector<Point3d> points = get_points(99000, 1);



  
  std::vector<Point3d> points;
  double m = 15.0;
  for(int i = 0; i < 30; i++) {
    for (int j = 0; j < 30; j++) {
      for(int k = 0; k < 30; k++) {
        double di = (double) i;
        double dj = (double) j;
        double dk = (double) k;
        Point3d p = {di-m, dj-m, dk-m};
        points.push_back(p);
      }
    }
  }
  
  
  
  

  Delaunay* del = generateDelaunay(points);
  std::vector<uint32_t> ttd = {0}; 
  draw_tetras_and_point(del, del->points[19]);
  return del;
}

void generateDelaunayTestMultiple() {
  for (int i = 0; i < 10000; i++) {
    std::cout << "Test : " << i << "\n";
    std::vector<Point3d> points = get_points(100, i);
    Delaunay* del = generateDelaunay(points);
  }
}



int do_delaunay() {

/*Point3d p0 = {1.0,0.0,0.0};
  Point3d p1 = {0.0,1.0,0.0};
  Point3d p2 = {0.0,0.0,1.0};
  Point3d p3 = {-1.0,0.0,0.0};
  Point3d p4 = {2.0,0.0,0.0};

  std::cout << "in shpere test: " << in_sphere_3d_SOS(p0,p1,p2,p3, p4) << "\n"; 
  std::cout << "in shpere test: " << in_sphere_3d_SOS(p1,p0,p2,p3, p4) << "\n"; 
  return 0;
  */
  

  try {
    app.initialize();
  }
  catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  generateDelaunayTest();
  return EXIT_SUCCESS;
  
}


