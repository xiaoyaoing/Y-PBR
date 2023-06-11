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

find_package(miniz REQUIRED)
find_package(Embree 3.0 REQUIRED)



## complie for Zero
file(GLOB_RECURSE YPBR_SOURCES ${YPBR_SOURCES_DIR}/*.cpp  lib/lodepng/lodepng.cpp lib/stbi/stb_image.c
        lib/skylight/ArHosekSkyModel.cpp)
link_libraries(spdlog::spdlog )
link_libraries(${EMBREE_LIBRARIES})


add_executable(Y_PBR test/main.cpp ${YPBR_SOURCES})
target_include_directories(Y_PBR PUBLIC ${PROJECT_SOURCE_DIR}/src )
target_link_libraries(Y_PBR PRIVATE tinyexr)


add_executable(sampleTest test/samplerTest.cpp  ${YPBR_SOURCES})
target_include_directories(sampleTest PUBLIC ${PROJECT_SOURCE_DIR}/src )
target_link_libraries(sampleTest PRIVATE tinyexr)