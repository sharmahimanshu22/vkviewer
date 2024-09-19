cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(vkviewer DESCRIPTION "VKViewer" LANGUAGES CXX)
add_executable(${PROJECT_NAME})

set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

#####################


# specify the C++ standard
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)


# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
	# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	endif()
	# reserve enough stack size
	target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
	# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
else()
	# set standard optimization level
	target_compile_options(${PROJECT_NAME} PUBLIC -O2)
	# reserve enough stack size	
	target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
  		target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	endif()
endif()



target_include_directories(${PROJECT_NAME} PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Vulkan REQUIRED)

add_subdirectory(deps/glfw)
#add_subdirectory(deps/geogram)
add_subdirectory(deps/glm)
add_subdirectory(deps/tinyobjloader)
# The one below is a single file external header library
include_directories(deps/stbimage) 

# For IMGUI, the recommended method is to just include imgui .cpp files in your target like your own source files
set(IMGUI_DIR deps/imgui)
include_directories(${IMGUI_DIR} ${IMGUI_DIR}/backends)


add_subdirectory(deps/openstl)


target_include_directories(${PROJECT_NAME} PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/deps/IndirectPredicates/include)

target_sources(${PROJECT_NAME} PUBLIC ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp ${IMGUI_DIR}/backends/imgui_impl_vulkan.cpp ${IMGUI_DIR}/imgui.cpp ${IMGUI_DIR}/imgui_draw.cpp ${IMGUI_DIR}/imgui_demo.cpp ${IMGUI_DIR}/imgui_tables.cpp ${IMGUI_DIR}/imgui_widgets.cpp)



# Final compilation and linking
target_sources(${PROJECT_NAME} PRIVATE main.cpp src/objLoader.cpp src/Predicates_psm.cpp src/widgets/widgets.cpp src/widgets/splinewidget.cpp src/voronoi3d.cpp src/utils.cpp src/vkviewer.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC glfw ${GLFW_LIBRARIES} Vulkan::Vulkan glm::glm tinyobjloader openstl_core) 
