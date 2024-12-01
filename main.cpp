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

#include "csg.h"

int main() {

  //do_csg();
  do_delaunay();
    
  return EXIT_SUCCESS;
}


