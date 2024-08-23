#include "objLoader.h"
#include "utils.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>


#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <iostream>

#include <openstl/core/stl.h>
#include "voronoi3d.h"

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



  DataForGPU loadSTL(std::string filename) {

    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
      std::cerr << "Error: Unable to open file '" << filename << "'" << std::endl;
    }
    
    // Deserialize the triangles in either binary or ASCII format
    std::vector<openstl::Triangle> triangles = openstl::deserializeStl(file);
    file.close();

    std::unordered_map<Vertex, uint32_t> uniqueVertices{};
    DataForGPU dataForGPU{};

    const auto& [vertices, faces] = convertToVerticesAndFaces(triangles);

    //std::cout << vertices.size() << "  vertices size\n";
    //std::cout << faces.size() << "  faces size\n";

    int j = 0;
    double xmax = -10000.0;
    double ymax = -10000.0;
    double zmax = -10000.0;
    double xmin = 10000.0;
    double ymin = 10000.0;
    double zmin = 10000.0;
    for(const auto& f : faces) {
      for(const auto& i : f) {

	Vertex vertex{};

	if (vertices[i].x > xmax) {
	  xmax = vertices[i].x;
	}
	if (vertices[i].y > ymax) {
	  ymax = vertices[i].y;
	}
	if (vertices[i].z > zmax) {
	  zmax = vertices[i].z;
	}
	if(vertices[i].x < xmin) {
	  xmin = vertices[i].x;
	}
	if(vertices[i].y < ymin) {
	  ymin = vertices[i].y;
	}
	if(vertices[i].z < zmin) {
	  zmin = vertices[i].z;
	}

	// Assuming that the mesh is triangular.
	vertex.pos = {vertices[i].x, vertices[i].y, vertices[i].z};
	vertex.color = {1.0f, 0.0f, 0.0f};
	vertex.texCoord = {0.0f,0.0f};
	
	if (uniqueVertices.count(vertex) == 0) {
	  uniqueVertices[vertex] = static_cast<uint32_t>(dataForGPU.vertices.size());
	  dataForGPU.vertices.push_back(vertex);
	}
	
	dataForGPU.indices.push_back(uniqueVertices[vertex]);
      }
    }
    std::cout << "maxmin\n" << xmax << " " << ymax << " " << zmax << "\n";
    std::cout << "maxmin\n" << xmin << " " << ymin << " " << zmin << "\n";

    std::cout << dataForGPU.indices.size() << " " << dataForGPU.vertices.size() << "  dataGPU\n";

    
    return dataForGPU;
  }



  DataForGPU loadDelaunay() {

    Delaunay del = generateDelaunayTest();
        
    std::unordered_map<Vertex, uint32_t> uniqueVertices{};
    DataForGPU dataForGPU{};
    
    Vertex vertex{};
    glm::dvec3 pt;
    bool boundary = false;
    int j = 0;
    for(const glm::uvec4& tet : del.tetrahedra) {
      
      uint32_t triangleIdxList[12]  = {tet[1],tet[2],tet[3],tet[0], tet[3], tet[2], tet[3],tet[0],tet[1],tet[1], tet[0], tet[2]};
      for(int j = 0; j < 12; j++) {
	if (triangleIdxList[j] < 8) {
	  boundary = true;
	  break;
	}
      }

      if (boundary) {
	boundary = false;
	continue;
      }

      for(int i = 0 ; i < 12; i++) {
	std::cout << triangleIdxList[i] << "," ;
	pt = del.points[triangleIdxList[i]];
	vertex.pos = {pt.x , pt.y, pt.z};
	vertex.color = {1.0f, 0.0f, 0.0f};
	vertex.texCoord = {0.0f,0.0f};
	
	if (uniqueVertices.count(vertex) == 0) {
	  uniqueVertices[vertex] = static_cast<uint32_t>(dataForGPU.vertices.size());
	  dataForGPU.vertices.push_back(vertex);
	}
	dataForGPU.indices.push_back(uniqueVertices[vertex]);
      }
      std::cout << "\n";
    }

    std::cout << dataForGPU.indices.size() << " " << dataForGPU.vertices.size() << "  dataGPU\n";
    return dataForGPU;
  }

}
