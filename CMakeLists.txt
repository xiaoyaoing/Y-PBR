
cmake_minimum_required(VERSION 3.18)
project(Y-PBR 	LANGUAGES CXX CUDA C
        )
# 设置C++标准为C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)









if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()



if(MSVC)
set(CUDAToolkit_ROOT "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.1")
set(CMAKE_CUDA_COMPILER "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1/bin/nvcc.exe")
# Adapted from the CMake source code at https://github.com/Kitware/CMake/blob/master/Modules/FindCUDA/select_compute_arch.cmake
# Simplified to return a semicolon-separated list of the compute capabilities of installed devices
function(TCNN_AUTODETECT_CUDA_ARCHITECTURES OUT_VARIABLE)
    if(NOT TCNN_AUTODETECT_CUDA_ARCHITECTURES_OUTPUT)
        if(CMAKE_CUDA_COMPILER_LOADED) # CUDA as a language
            set(file "${PROJECT_BINARY_DIR}/detect_tcnn_cuda_architectures.cu")
        else()
            set(file "${PROJECT_BINARY_DIR}/detect_tcnn_cuda_architectures.cpp")
        endif()

        file(WRITE ${file} ""
                "#include <cuda_runtime.h>\n"
                "#include <cstdio>\n"
                "int main() {\n"
                "	int count = 0;\n"
                "	if (cudaSuccess != cudaGetDeviceCount(&count)) return -1;\n"
                "	if (count == 0) return -1;\n"
                "	for (int device = 0; device < count; ++device) {\n"
                "		cudaDeviceProp prop;\n"
                "		if (cudaSuccess == cudaGetDeviceProperties(&prop, device)) {\n"
                "			std::printf(\"%d%d\", prop.major, prop.minor);\n"
                "			if (device < count - 1) std::printf(\";\");\n"
                "		}\n"
                "	}\n"
                "	return 0;\n"
                "}\n"
                )

        if(CMAKE_CUDA_COMPILER_LOADED) # CUDA as a language
            try_run(run_result compile_result ${PROJECT_BINARY_DIR} ${file} RUN_OUTPUT_VARIABLE compute_capabilities)
        else()
            try_run(
                    run_result compile_result ${PROJECT_BINARY_DIR} ${file}
                    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CUDA_INCLUDE_DIRS}"
                    LINK_LIBRARIES ${CUDA_LIBRARIES}
                    RUN_OUTPUT_VARIABLE compute_capabilities
            )
        endif()

        if(run_result EQUAL 0)
            # If the user has multiple GPUs with the same compute capability installed, list that capability only once.
            list(REMOVE_DUPLICATES compute_capabilities)
            set(TCNN_AUTODETECT_CUDA_ARCHITECTURES_OUTPUT ${compute_capabilities} CACHE INTERNAL "Returned GPU architectures from detect_gpus tool" FORCE)
        endif()
    endif()

    if(NOT TCNN_AUTODETECT_CUDA_ARCHITECTURES_OUTPUT)
        message(STATUS "Automatic GPU detection failed. Building for Turing and Ampere as a best guess.")
        set(${OUT_VARIABLE} "75;86" PARENT_SCOPE)
    else()
        set(${OUT_VARIABLE} ${TCNN_AUTODETECT_CUDA_ARCHITECTURES_OUTPUT} PARENT_SCOPE)
    endif()
endfunction()

message(STATUS "Obtained CUDA architectures automatically from installed GPUs")
TCNN_AUTODETECT_CUDA_ARCHITECTURES(CMAKE_CUDA_ARCHITECTURES)


endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(YPBR_DEPENDENCIES_DIR ${PROJECT_SOURCE_DIR}/lib)
set(YPBR_SOURCES_DIR ${PROJECT_SOURCE_DIR}/src)


#include_directories(${PROJECT_SOURCE_DIR}/lib/json/)
include_directories(${PROJECT_SOURCE_DIR}/lib/loadpng/)
include_directories(${PROJECT_SOURCE_DIR}/lib/stbi/)


include_directories(${YPBR_SOURCE_DIR})



## complie for dependencies
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF)
add_subdirectory(${YPBR_DEPENDENCIES_DIR})





file(GLOB_RECURSE YPBR_SOURCES ${YPBR_SOURCES_DIR}/*.cpp  lib/lodepng/lodepng.cpp
        lib/skylight/ArHosekSkyModel.cpp)
include_directories(${PROJECT_SOURCE_DIR}/src)




add_executable(Y_PBR test/main.cpp ${YPBR_SOURCES})
add_dependencies(Y_PBR lib)
target_link_libraries(Y_PBR PUBLIC lib)
target_include_directories(Y_PBR PRIVATE ${PROJECT_SOURCE_DIR}/lib/json/include)

if(MSVC)
    message("Using MSVC compiler")
else()
    message("Not using MSVC compiler")
endif()

if(MSVC )
add_subdirectory(tiny-cuda-nn)

#target_link_libraries(Y_PBR PUBLIC E:/code/Y-PBR/lib/embree/windows/lib/embree3.lib)

add_executable(hair-nn test/hair-accleerate.cu ${YPBR_SOURCES}  test/HairIntegrator.cpp test/HairIntegrator.h src/Bsdfs/PrecomputeALobe.cpp src/Bsdfs/PrecomputeALobe.h)
target_compile_features(hair-nn PRIVATE cxx_std_17)

add_dependencies(hair-nn lib)
target_link_libraries(hair-nn PUBLIC lib ${CUDA_LIBRARIES} tiny-cuda-nn)
target_include_directories(hair-nn PRIVATE ${PROJECT_SOURCE_DIR}/lib/json/)

endif()


add_executable(sampleTest test/samplerTest.cpp  ${YPBR_SOURCES})
target_link_libraries(sampleTest PRIVATE tinyexr)

if (WIN32)
    install(TARGETS ${executables} DESTINATION .)
else()
    install(TARGETS ${executables} DESTINATION bin)
endif()
