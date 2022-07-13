cmake_minimum_required(VERSION 3.10)
project(Y-PBR)


set (CMAKE_C_COMPILER "/usr/local/bin/g++-11")
set (CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)




include_directories(${PROJECT_SOURCE_DIR}/lib/glm/)
include_directories(${PROJECT_SOURCE_DIR}/lib/nlohmann/)
include_directories(${PROJECT_SOURCE_DIR}/lib/spdlog/)


file(GLOB_RECURSE _source_list ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.hpp )




foreach(_source IN ITEMS ${_source_list})
    get_filename_component(_source_path "${_source}" PATH)
    string(REPLACE "${CMAKE_SOURCE_DIR}" "" _group_path "${_source_path}")
    string(REPLACE "/" "\\" _group_path "${_group_path}")
    source_group("${_group_path}" FILES "${_source}")
endforeach()

message(${_source_list})

add_executable(${PROJECT_NAME} ${_source_list})
