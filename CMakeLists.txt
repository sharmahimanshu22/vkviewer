cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(vkviewer DESCRIPTION "VKViewer" LANGUAGES CXX)


add_executable(${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Vulkan REQUIRED)

add_subdirectory(deps/glfw)
add_subdirectory(deps/geogram)
add_subdirectory(deps/glm)
add_subdirectory(deps/tinyobjloader)
# The one below is a single file external header library
include_directories(deps/stbimage) 

# For IMGUI, the recommended method is to just include imgui .cpp files in your target like your own source files
set(IMGUI_DIR deps/imgui)
include_directories(${IMGUI_DIR} ${IMGUI_DIR}/backends)
target_sources(${PROJECT_NAME} PUBLIC ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp ${IMGUI_DIR}/backends/imgui_impl_vulkan.cpp ${IMGUI_DIR}/imgui.cpp ${IMGUI_DIR}/imgui_draw.cpp ${IMGUI_DIR}/imgui_demo.cpp ${IMGUI_DIR}/imgui_tables.cpp ${IMGUI_DIR}/imgui_widgets.cpp)


add_subdirectory(deps/openstl)

# Final compilation and linking
target_sources(${PROJECT_NAME} PRIVATE main.cpp src/objLoader.cpp )
target_link_libraries(${PROJECT_NAME} PUBLIC glfw ${GLFW_LIBRARIES} Vulkan::Vulkan glm::glm tinyobjloader openstl_core geogram) 
