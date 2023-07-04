#add_subdirectory(./glm)
#add_subdirectory(./catch2)
#add_subdirectory(./tinyexr)
#add_subdirectory(./spdlog)
#add_subdirectory(./embree)

set(TARGET_NAME lib)
set(3DPARTY_FOLDER lib)

function(set_ide_folder target_name folder_name)
    if (TARGET ${target_name})
        set_target_properties(${target_name} PROPERTIES FOLDER ${folder_name})
    endif()
endfunction()


add_library(${TARGET_NAME} SHARED lib.cpp)



if (NOT TARGET embree)
    message(STATUS "============start config embree============")
    set(EMBREE_FOLDER ${3DPARTY_FOLDER}/embree)
    if(EMBREE_USE_TBB)
        set(EMBREE_TASKING_SYSTEM TBB)
    else()
        if(WIN32)
            set(EMBREE_TASKING_SYSTEM PPL)
        else()
            set(EMBREE_TASKING_SYSTEM INTERNAL)
        endif()
    endif()

    set(EMBREE_TUTORIALS OFF CACHE INTERNAL "")

    # add_subdirectory(embree)
    # target_link_libraries(${TARGET_NAME} PUBLIC embree)
    # set_target_properties(embree PROPERTIES DEBUG_POSTFIX "_d")
    # set_target_properties(embree PROPERTIES RELEASE_POSTFIX "_r")

    set(EMBREE_VERSION "3.13.5")
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/embree)
    message(STATUS "Downloading embree ${EMBREE_VERSION}...")
    include(FetchContent)
    if (WIN32)
        set(EMBREE_URL_SUFFIX ".x64.vc14.windows.zip")
    elseif (APPLE)
        set(EMBREE_URL_SUFFIX ".x86_64.macosx.zip")
    else ()# Linux
        set(EMBREE_URL_SUFFIX ".x86_64.linux.tar.gz")
    endif ()

    set(EMBREE_DOWNLOAD_URL "https://github.com/embree/embree/releases/download/v${EMBREE_VERSION}/embree-${EMBREE_VERSION}${EMBREE_URL_SUFFIX}")
    message(STATUS ${embree_library})
    FetchContent_Declare(
            embree_binary
            URL ${EMBREE_DOWNLOAD_URL}
            GIT_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(embree_binary)
    FetchContent_GetProperties(embree_binary)
    # patch embree
    message(STATUS "Patching embree-config.cmake...")
    file(READ "${embree_binary_SOURCE_DIR}/lib/cmake/embree-${EMBREE_VERSION}/embree-config.cmake" EMBREE_CONFIG)
    string(REPLACE "FIND_PACKAGE(TBB REQUIRED \${EMBREE_TBB_COMPONENT})" "" EMBREE_CONFIG "${EMBREE_CONFIG}")
    file(WRITE "${embree_binary_SOURCE_DIR}/lib/cmake/embree-${EMBREE_VERSION}/embree-config.cmake" "${EMBREE_CONFIG}")
    message(STATUS "Downloaded embree: ${embree_binary_SOURCE_DIR}")
    set(embree_ROOT ${embree_binary_SOURCE_DIR})
    message("embree_root" ${embree_ROOT})

    #    list(CMAKE_PREFIX_PATH APPEND ${embree_ROOT} )
    find_package(embree 3 CONFIG REQUIRED)

    #    add_custom_command(TARGET lib-copy
    #            DEPENDS $<TARGET_FILE:embree>
    #            COMMAND ${CMAKE_COMMAND} -E copy_if_different
    #            $<TARGET_FILE:embree> $<TARGET_FILE_DIR:${TARGET_NAME}>)
    #    add_custom_command(TARGET lib-copy
    #            DEPENDS $<TARGET_LINKER_FILE:embree>
    #            COMMAND ${CMAKE_COMMAND} -E copy_if_different
    #            $<TARGET_LINKER_FILE:embree> $<TARGET_FILE_DIR:${TARGET_NAME}>)
    target_include_directories(embree INTERFACE ${embree_binary_SOURCE_DIR}/include)
    target_link_libraries(${TARGET_NAME} INTERFACE embree)
endif()

if(NOT TARGET glm)
    message(STATUS "============start config glm============")
    add_subdirectory(glm)
    target_link_libraries(${TARGET_NAME} INTERFACE glm )
    target_include_directories(glm INTERFACE ${CMAKE_CURRENT_LIST_DIR}/glm)
endif ()


target_link_libraries(${TARGET_NAME} INTERFACE libatomic)

if(NOT TARGET spdlog)
    #spdlog is header-only library.So just include
    #add_subdirectory(spdlog)
   # target_link_libraries(${TARGET_NAME} INTERFACE spdlog )
    #target_include_directories(spdlog PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spdlog/include)
    target_include_directories(${TARGET_NAME} INTERFACE spdlog)
endif ()



if(NOT TARGET tinyexr)
    message(STATUS "============start config glm============")
    add_subdirectory(tinyexr)
    target_link_libraries(${TARGET_NAME} INTERFACE tinyexr)
    target_include_directories(tinyexr INTERFACE ${CMAKE_CURRENT_LIST_DIR}/tinyexr)
endif ()


