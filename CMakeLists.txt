cmake_minimum_required(VERSION 3.10)
project(Y-PBR)


set( CMAKE_BUILD_TYPE Debug )


#set (CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
#set (CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(YPBR_DEPENDENCIES_DIR ${PROJECT_SOURCE_DIR}/lib)
set(YPBR_SOURCES_DIR ${PROJECT_SOURCE_DIR}/src)


include_directories(${PROJECT_SOURCE_DIR}/lib/nlohmann/)
include_directories(${PROJECT_SOURCE_DIR}/lib/loadpng/)
include_directories(${PROJECT_SOURCE_DIR}/lib/stbi/)



include_directories(${YPBR_DEPENDENCIES_DIR})
include_directories(${YPBR_SOURCE_DIR})
include_directories($<BUILD_INTERFACE:${EMBREE_INCLUDE_DIRS}>)



## complie for dependencies
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF)
add_subdirectory(${YPBR_DEPENDENCIES_DIR})


if(NOT TARGET spdlog)
    find_package(spdlog REQUIRED)
endif()

find_package(Embree 3.0 REQUIRED)

## complie for Zero
file(GLOB_RECURSE YPBR_SOURCES ${YPBR_SOURCES_DIR}/*.cpp src/Common/Parallel.cpp lib/lodepng/lodepng.cpp lib/stbi/stb_image.c
        lib/skylight/ArHosekSkyModel.cpp)
link_libraries(spdlog::spdlog )
link_libraries(${EMBREE_LIBRARIES})

add_executable(Y_PBR src/main.cpp ${YPBR_SOURCES} src/Integrator/BDPT/PathVertex.cpp src/Integrator/BDPT/PathVertex.h src/Integrator/BDPT/Records.h src/Integrator/BDPT/LightPath.cpp src/Integrator/BDPT/LightPath.h src/Integrator/BDPT/CameraPath.cpp src/Integrator/BDPT/CameraPath.h src/Integrator/BDPT/BdptTracer.cpp src/Integrator/BDPT/BdptTracer.h src/Bsdfs/RoughHair.h src/Bsdfs/RoughHair.cpp src/Integrator/TraceHelper.cpp src/Integrator/TraceHelper.h src/SampleRecords/records.h src/SampleRecords/PositionAndDirectionSample.h src/Camera/ImagePramid.cpp src/Camera/ImagePramid.h)
target_include_directories(Y_PBR PUBLIC ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/lib)



add_subdirectory(./test)
#
#
#
#
#
#if(NOT TARGET spdlog)
#    find_package(spdlog REQUIRED)
#endif()
#
#
#include_directories(${PROJECT_SOURCE_DIR}/lib/glm/)
#include_directories(${PROJECT_SOURCE_DIR}/lib/nlohmann/)
#
#
#enable_testing()
#add_subdirectory(lib/gtest)
#include_directories(lib/gtest)
#add_tests(tests)
#
#
#file(GLOB_RECURSE _source_list ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.hpp )
#
#
#foreach(_source IN ITEMS ${_source_list})
#    get_filename_component(_source_path "${_source}" PATH)
#    string(REPLACE "${CMAKE_SOURCE_DIR}" "" _group_path "${_source_path}")
#    string(REPLACE "/" "\\" _group_path "${_group_path}")
#    source_group("${_group_path}" FILES "${_source}")
#endforeach()
##
#message(${_source_list})
#
#
#add_executable(${PROJECT_NAME} ${_source_list} src/Primitives/Cube.hpp  src/Primitives/TriangleMesh.hpp src/Primitives/TriangleMesh.cpp )
#target_link_libraries(${PROJECT_NAME}  spdlog::spdlog)

