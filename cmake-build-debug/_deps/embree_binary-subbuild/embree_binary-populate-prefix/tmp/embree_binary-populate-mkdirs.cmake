# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-src"
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-build"
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix"
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix/tmp"
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix/src/embree_binary-populate-stamp"
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix/src"
  "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix/src/embree_binary-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix/src/embree_binary-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "E:/code/Y-PBR/cmake-build-debug/_deps/embree_binary-subbuild/embree_binary-populate-prefix/src/embree_binary-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
