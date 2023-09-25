# add_subdirectory(./glm)
# add_subdirectory(./catch2)
# add_subdirectory(./tinyexr)
# add_subdirectory(./spdlog)
# add_subdirectory(./embree)
cmake_minimum_required(VERSION 3.14)

set(TARGET_NAME lib)
set(3DPARTY_FOLDER lib)
set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Build Embree with support for ISPC applications." FORCE)

function(set_ide_folder target_name folder_name)
    if(TARGET ${target_name})
        set_target_properties(${target_name} PROPERTIES FOLDER ${folder_name})
    endif()
endfunction()

add_library(${TARGET_NAME} STATIC lib.cpp)

if(WIN32)
    set(SYSTEM_TYPE "windows")
elseif(APPLE)
    set(SYSTEM_TYPE "macos")
else() # Linux
    message(STATUS "linuxe not support")
endif()

set(CMAKE_PREFIX_PATH APPEND "${CMAKE_CURRENT_SOURCE_DIR}/embree/${SYSTEM_TYPE}")
set(EMBREE_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/embree/${SYSTEM_TYPE}")
message(STATUS ${EMBREE_LIB_DIR})
find_package(embree 3 REQUIRED)

target_link_libraries(${TARGET_NAME} INTERFACE embree)
target_include_directories(embree INTERFACE ${CMAKE_CURRENT_LIST_DIR}/embree/include)

message(STATUS ${EMBREE_LIB_DIR}/bin/embree3.dll)
message(STATUS ${CMAKE_BINARY_DIR}/tbb12.dll)

if(WIN32)
    file(COPY ${EMBREE_LIB_DIR}/bin/embree3.dll DESTINATION ${CMAKE_BINARY_DIR})
    file(COPY ${EMBREE_LIB_DIR}/bin/tbb12.dll DESTINATION ${CMAKE_BINARY_DIR})
endif()


if(NOT TARGET glm)
    message(STATUS "============start config glm============")
    add_subdirectory(glm)
    target_include_directories(glm INTERFACE ${CMAKE_CURRENT_LIST_DIR}/glm)
    target_link_libraries(${TARGET_NAME} INTERFACE glm)
endif()

if(NOT MSVC AND WIN32)
    target_link_libraries(${TARGET_NAME} INTERFACE libatomic)
endif()

if(NOT TARGET spdlog)
    # spdlog is header-only library.So just include
    # add_subdirectory(spdlog)
    # target_link_libraries(${TARGET_NAME} INTERFACE spdlog )
    # target_include_directories(spdlog PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spdlog/include)
    target_include_directories(${TARGET_NAME} INTERFACE spdlog)
endif()

if(NOT TARGET tinyexr)
    message(STATUS "============start config glm============")
    add_subdirectory(tinyexr)
    target_link_libraries(${TARGET_NAME} INTERFACE tinyexr)
    target_include_directories(tinyexr INTERFACE ${CMAKE_CURRENT_LIST_DIR}/tinyexr)
endif()

if(NOT TARGET stb)
    message(STATUS "============start config stb============")
    add_subdirectory(stb)
    target_link_libraries(${TARGET_NAME} PUBLIC stb)
    target_include_directories(stb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/stb)
endif()

# add_subdirectory(json)
# target_link_libraries( ${TARGET_NAME} PUBLIC nlohmann_json)
#file(COPY ${EMBREE_LIB_DIR} DESTINATION ${CMAKE_BINARY_DIR})

