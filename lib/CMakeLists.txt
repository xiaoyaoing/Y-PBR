#add_subdirectory(./glm)
#add_subdirectory(./catch2)
#add_subdirectory(./tinyexr)
#add_subdirectory(./spdlog)
#add_subdirectory(./embree)
cmake_minimum_required(VERSION 3.14)

set(TARGET_NAME lib)
set(3DPARTY_FOLDER lib)
set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Build Embree with support for ISPC applications." FORCE)

function(set_ide_folder target_name folder_name)
    if (TARGET ${target_name})
        set_target_properties(${target_name} PROPERTIES FOLDER ${folder_name})
    endif ()
endfunction()


add_library(${TARGET_NAME} STATIC lib.cpp)

if (WIN32)
    set(SYSTEM_TYPE "windows")
elseif (APPLE)
    set(SYSTEM_TYPE "macos")
else ()# Linux
    message(STATUS "linuxe not support")
endif ()


set(CMAKE_PREFIX_PATH APPEND "${CMAKE_CURRENT_SOURCE_DIR}/embree/${SYSTEM_TYPE}")
set(EMBREE_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/embree/${SYSTEM_TYPE}")
message(STATUS ${EMBREE_LIB_DIR} )
find_package(embree 3  REQUIRED)

target_link_libraries(${TARGET_NAME} INTERFACE embree )
target_include_directories(embree INTERFACE ${CMAKE_CURRENT_LIST_DIR}/embree/include)

if(WIN32)
add_custom_command(TARGET lib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${EMBREE_LIB_DIR}/bin/embree3.dll $<TARGET_FILE_DIR:${TARGET_NAME}>/embree3.dll

        COMMAND ${CMAKE_COMMAND} -E copy
        ${EMBREE_LIB_DIR}/bin/tbb12.dll $<TARGET_FILE_DIR:${TARGET_NAME}>/tbb12.dll)

endif()

if (NOT TARGET glm)
    message(STATUS "============start config glm============")
    add_subdirectory(glm)
    target_link_libraries(${TARGET_NAME} INTERFACE glm)
    target_include_directories(glm INTERFACE ${CMAKE_CURRENT_LIST_DIR}/glm)
endif ()

if(NOT MSVC AND WIN32)
    target_link_libraries(${TARGET_NAME} INTERFACE libatomic)
endif()


if (NOT TARGET spdlog)
    #spdlog is header-only library.So just include
    #add_subdirectory(spdlog)
    # target_link_libraries(${TARGET_NAME} INTERFACE spdlog )
    #target_include_directories(spdlog PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spdlog/include)
    target_include_directories(${TARGET_NAME} INTERFACE spdlog)
endif ()


if (NOT TARGET tinyexr)
    message(STATUS "============start config glm============")
    add_subdirectory(tinyexr)
    target_link_libraries(${TARGET_NAME} INTERFACE tinyexr)
    target_include_directories(tinyexr INTERFACE ${CMAKE_CURRENT_LIST_DIR}/tinyexr)
endif ()

#add_subdirectory(json)
#target_link_libraries( ${TARGET_NAME} PUBLIC nlohmann_json)


