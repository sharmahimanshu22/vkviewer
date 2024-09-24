#include "vkviewer.h"
#include "objLoader.h"
#include <iostream>
#include "voronoi3d.h"
#include <cstdlib>

#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/ext/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale
#include <glm/ext/matrix_clip_space.hpp> // glm::perspective
#include <glm/ext/scalar_constants.hpp> // glm::pi
#include <iostream>

std::ostream &operator<<(std::ostream &os, glm::mat4 const &m) {
 return os << m[0][0] << " " << m[1][0] << " " << m[2][0] << " " << m[3][0] << "\n" << m[0][1] << " " << m[1][1] << " " << m[2][1] << " " << m[3][1] << "\n"  << m[0][2] << " " << m[1][2] << " " << m[2][2] << " " << m[3][2] << "\n" << m[0][3] << " " << m[1][3] << " " << m[2][3] << " " << m[3][3] ;
}

glm::mat4 camera(float Translate, glm::vec2 const& Rotate)
{
	glm::mat4 Projection = glm::perspective(glm::pi<float>() * 0.25f, 4.0f / 3.0f, 1.0f, 100.f);
	//std::cout << Projection << "\n";
	glm::mat4 View = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -Translate));
	View = glm::rotate(View, Rotate.y, glm::vec3(-1.0f, 0.0f, 0.0f));
	View = glm::rotate(View, Rotate.x, glm::vec3(0.0f, 1.0f, 0.0f));
	glm::mat4 Model = glm::scale(glm::mat4(1.0f), glm::vec3(0.5f));
	return Projection * View * Model;
}


int main() {

  //camera(0.1, glm::vec2(1.0,1.0));
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


