cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(myapp DESCRIPTION "VKViewer" LANGUAGES CXX)

#enable_language(CUDA)
#if(NOT DEFINED CMAKE_CUDA_STANDARD)
#    set(CMAKE_CUDA_STANDARD 17)
#    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
#endif()

#if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
#    set(CMAKE_CUDA_ARCHITECTURES 61)
#endif()

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE main.cpp src/objLoader.cpp)

#set_target_properties(myapp PROPERTIES
#                            CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(${PROJECT_NAME}
        PUBLIC
                ${CMAKE_CURRENT_SOURCE_DIR}/include
)

#include_directories(/Users/sharmh15/VulkanSDK/1.3.275.0/macOS/include)

#include_directories(include)

find_package(Vulkan REQUIRED)

add_subdirectory(deps/glfw)
add_subdirectory(deps/glm)
add_subdirectory(deps/tinyobjloader)
# The one below is a single file external header library
include_directories(deps/stbimage) 

target_link_libraries(${PROJECT_NAME} PUBLIC glfw ${GLFW_LIBRARIES} Vulkan::Vulkan glm::glm tinyobjloader) 



#add_subdirectory(deps/tinyobjloader)
#target_link_libraries(${PROJECT_NAME} PUBLIC tinyobjloader)
