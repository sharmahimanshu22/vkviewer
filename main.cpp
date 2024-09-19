#include "vkviewer.h"
#include "objLoader.h"
#include <iostream>
#include "voronoi3d.h"
#include <cstdlib>

int main() {

  do_delaunay();
  
  //double p0[2] = {0.0,0.0};
  //double p1[2] = {1.0,0.0};
  //double p2[2] = {1.0,-1.0};
  
  //std::cout << GEO::PCK::orient_2d(p0,p1,p2) << "\n";

  //generateDelaunayTest();
  //vkview::loadSTL("/Users/sharmh15/projectsPersonal/Thingi10K/raw_meshes/242237.stl");
  
  //return 0;

  /*
  HelloTriangleApplication app;

  vkview::DataForGPU dfg = vkview::loadDelaunay();
  
  try {
    app.initialize();
    app.loadScene(dfg);
    app.run();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  */
  
  return EXIT_SUCCESS;
}


