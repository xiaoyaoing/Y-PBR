cmake_minimum_required(VERSION 3.19)
project(Y-PBR)


if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()



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





file(GLOB_RECURSE YPBR_SOURCES ${YPBR_SOURCES_DIR}/*.cpp  lib/lodepng/lodepng.cpp lib/stbi/stb_image.c
        lib/skylight/ArHosekSkyModel.cpp)
include_directories(${PROJECT_SOURCE_DIR}/src)




add_executable(Y_PBR test/main.cpp ${YPBR_SOURCES})



add_dependencies(Y_PBR lib)
target_link_libraries(Y_PBR PUBLIC lib)

#target_link_libraries(Y_PBR PUBLIC E:/code/Y-PBR/lib/embree/windows/lib/embree3.lib)



#
#set(TARGET_NAME Y_PBR)
#add_custom_command(TARGET Y_PBR POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_CURRENT_SOURCE_DIR}/lib/embree/bin/embree3.dll $<TARGET_FILE_DIR:${TARGET_NAME}>/embree3.dll
#
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_CURRENT_SOURCE_DIR}/lib/embree/bin/tbb12.dll $<TARGET_FILE_DIR:${TARGET_NAME}>/tbb12.dll
#        )
##



add_executable(sampleTest test/samplerTest.cpp  ${YPBR_SOURCES})
target_link_libraries(sampleTest PRIVATE tinyexr)

if (WIN32)
    install(TARGETS ${executables} DESTINATION .)
else()
    install(TARGETS ${executables} DESTINATION bin)
endif()
