#ifndef OBJLOADER
#define OBJLOADER

#include <vector>
#include "vertex.h"


namespace vkview {
  
  struct DataForGPU {
    std::vector<vkview::Vertex> vertices;
    std::vector<uint32_t> indices;
    
  };

  DataForGPU loadModel(const std::string model_path);

  DataForGPU loadSTL(const std::string filename);

  DataForGPU loadDelaunay();
}

#endif
