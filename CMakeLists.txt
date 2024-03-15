
cmake_minimum_required(VERSION 3.18)

message(STATUS "CMake version: ${CMAKE_VERSION}")


project(Y-PBR 	LANGUAGES CXX  C)
# 设置C++标准为C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()




set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(YPBR_DEPENDENCIES_DIR ${PROJECT_SOURCE_DIR}/lib)
set(YPBR_SOURCES_DIR ${PROJECT_SOURCE_DIR}/src)


include_directories(${PROJECT_SOURCE_DIR}/lib/nlohmann/)
include_directories(${PROJECT_SOURCE_DIR}/lib/loadpng/)
#include_directories(${PROJECT_SOURCE_DIR}/lib/stbi/)


include_directories(${YPBR_SOURCE_DIR})



## complie for dependencies
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF)
add_subdirectory(${YPBR_DEPENDENCIES_DIR})





file(GLOB_RECURSE YPBR_SOURCES ${YPBR_SOURCES_DIR}/*.cpp   ${YPBR_SOURCES_DIR}/*.h  lib/lodepng/lodepng.cpp ${YPBR_SOURCES_DIR}/*.hpp
        lib/skylight/ArHosekSkyModel.cpp)
include_directories(${PROJECT_SOURCE_DIR}/src)




add_executable(Y_PBR test/main.cpp ${YPBR_SOURCES})
add_dependencies(Y_PBR lib)
target_link_libraries(Y_PBR PUBLIC lib)
#target_include_directories(Y_PBR PRIVATE ${PROJECT_SOURCE_DIR}/lib/nlohmann/include)

if(MSVC)
    message("Using MSVC compiler")
else()
    message("Not using MSVC compiler")
endif()





#target_link_libraries(Y_PBR PUBLIC E:/code/Y-PBR/lib/embree/windows/lib/embree3.lib)

#target_include_directories(hair-nn PRIVATE ${PROJECT_SOURCE_DIR}/lib/nlohmann/)



add_executable(sampleTest test/samplerTest.cpp  ${YPBR_SOURCES})
target_link_libraries(sampleTest PRIVATE tinyexr)

if (WIN32)
    install(TARGETS ${executables} DESTINATION .)
else()
    install(TARGETS ${executables} DESTINATION bin)
endif()
