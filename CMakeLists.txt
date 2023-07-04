cmake_minimum_required(VERSION 3.19)
project(Y-PBR)



if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()

#set (CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
#set (CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(YPBR_DEPENDENCIES_DIR ${PROJECT_SOURCE_DIR}/lib)
set(YPBR_SOURCES_DIR ${PROJECT_SOURCE_DIR}/src)


include_directories(${PROJECT_SOURCE_DIR}/lib/nlohmann/)
include_directories(${PROJECT_SOURCE_DIR}/lib/loadpng/)
include_directories(${PROJECT_SOURCE_DIR}/lib/stbi/)


include_directories(${YPBR_SOURCE_DIR})



## complie for dependencies
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF)
add_subdirectory(${YPBR_DEPENDENCIES_DIR})

#
#if(NOT TARGET spdlog)
#    find_package(spdlog REQUIRED)
#endif()

#find_package(miniz REQUIRED)



file(GLOB_RECURSE YPBR_SOURCES ${YPBR_SOURCES_DIR}/*.cpp  lib/lodepng/lodepng.cpp lib/stbi/stb_image.c
        lib/skylight/ArHosekSkyModel.cpp)
#link_libraries(spdlog::spdlog )
include_directories(${PROJECT_SOURCE_DIR}/src)




add_executable(Y_PBR test/main.cpp ${YPBR_SOURCES})



add_dependencies(Y_PBR lib)
target_link_libraries(Y_PBR PUBLIC lib)

#target_include_directories(Y_PBR PUBLIC  )



add_executable(sampleTest test/samplerTest.cpp  ${YPBR_SOURCES})
#target_include_directories(sampleTest PUBLIC ${PROJECT_SOURCE_DIR}/src )
target_link_libraries(sampleTest PRIVATE tinyexr)


#config spdlog begin
#add_subdirectory(lib/spdlog)
#target_include_directories(Y_PBR PUBLIC lib/spdlog/include)
#target_link_libraries(Y_PBR INTERFACE spdlog)
#config spdlog end

#config spdlog begin
#add_subdirectory(lib/glm)
#target_include_directories(Y_PBR PUBLIC lib/glm)
#target_link_libraries(Y_PBR INTERFACE glm)
#config spdlog end