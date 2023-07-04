# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.25.2)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(embree_binary-populate NONE)



include(ExternalProject)
ExternalProject_Add(embree_binary-populate
                     "UPDATE_DISCONNECTED" "False" "DOWNLOAD_EXTRACT_TIMESTAMP" "FALSE" "URL" "https://github.com/embree/embree/releases/download/v3.13.5/embree-3.13.5.x64.vc14.windows.zip" "GIT_PROGRESS" "TRUE"
                    SOURCE_DIR          "E:/code/Y-PBR/cmake-build-release/_deps/embree_binary-src"
                    BINARY_DIR          "E:/code/Y-PBR/cmake-build-release/_deps/embree_binary-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


