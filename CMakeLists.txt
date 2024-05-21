cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(vkviewer DESCRIPTION "VKViewer" LANGUAGES CXX)


add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE main.cpp src/objLoader.cpp)

target_include_directories(${PROJECT_NAME}
        PUBLIC
                ${CMAKE_CURRENT_SOURCE_DIR}/include
)


find_package(Vulkan REQUIRED)

add_subdirectory(deps/glfw)
add_subdirectory(deps/glm)
add_subdirectory(deps/tinyobjloader)
# The one below is a single file external header library
include_directories(deps/stbimage) 

target_link_libraries(${PROJECT_NAME} PUBLIC glfw ${GLFW_LIBRARIES} Vulkan::Vulkan glm::glm tinyobjloader) 


