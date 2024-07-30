#include "objLoader.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>


#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <iostream>


namespace vkview {
  
  DataForGPU loadModel(const std::string model_path) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;
    DataForGPU dataForGPU{};
    
    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, model_path.c_str())) {
      throw std::runtime_error(warn + err);
    }
    
    
    std::unordered_map<Vertex, uint32_t> uniqueVertices{};
    
    for (const auto& shape : shapes) {
      for (const auto& index : shape.mesh.indices) {
	Vertex vertex{};
	
	// Assuming that the mesh is triangular.
	vertex.pos = {
	  attrib.vertices[3 * index.vertex_index + 0],
	  attrib.vertices[3 * index.vertex_index + 1],
	  attrib.vertices[3 * index.vertex_index + 2]
	};
	
	vertex.texCoord = {
	  attrib.texcoords[2 * index.texcoord_index + 0],
	  1.0f - attrib.texcoords[2 * index.texcoord_index + 1]    // Subtracting from 1.0f as texcoord and images are inverted in y index.
	};
	
	vertex.color = {1.0f, 1.0f, 1.0f};
	
	if (uniqueVertices.count(vertex) == 0) {
	  uniqueVertices[vertex] = static_cast<uint32_t>(dataForGPU.vertices.size());
	  dataForGPU.vertices.push_back(vertex);
	}
	
	dataForGPU.indices.push_back(uniqueVertices[vertex]);
      }
    }
    std::cout << dataForGPU.indices.size() << " should have worked\n";
    return dataForGPU;
  }


}
