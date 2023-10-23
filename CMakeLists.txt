 # ============================================================================
 #  This file is part of ASLI (A Simple Lattice Infiller)
 #  Copyright (C) KU Leuven, 2019-2022
 #
 #  ASLI is free software: you can redistribute it and/or modify it under the 
 #  terms of the GNU Affero General Public License as published by the Free 
 #  Software Foundation, either version 3 of the License, or (at your option) 
 #  any later version.
 #
 #  ASLI is distributed in the hope that it will be useful, but WITHOUT ANY 
 #  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 #  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for 
 #  more details.
 #
 #  You should have received a copy of the GNU Affero General Public License
 #  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 #  Please read the terms carefully and use this copy of ASLI only if you
 #  accept them.
 # ============================================================================

cmake_minimum_required(VERSION 3.13)

project(ASLI)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

include(GNUInstallDirs) # Variables pointing to standard install directories
include(ExternalProject)

# Release version and date
file(READ "include/version.h" ver)

string(REGEX MATCH "ASLI_VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ASLI_VERSION_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "ASLI_VERSION_MINOR ([0-9]*)" _ ${ver})
set(ASLI_VERSION_MINOR ${CMAKE_MATCH_1})
set(ASLI_VERSION "${ASLI_VERSION_MAJOR}.${ASLI_VERSION_MINOR}")

string(REGEX MATCH "ASLI_RELEASE_DAY ([0-9]*)" _ ${ver})
set(ASLI_RELEASE_DAY ${CMAKE_MATCH_1})
string(REGEX MATCH "ASLI_RELEASE_MONTH ([(A-z)|(a-z)]*)" _ ${ver})
set(ASLI_RELEASE_MONTH ${CMAKE_MATCH_1})
string(REGEX MATCH "ASLI_RELEASE_YEAR ([0-9]*)" _ ${ver})
set(ASLI_RELEASE_YEAR ${CMAKE_MATCH_1})
set(ASLI_RELEASE_DATE "${ASLI_RELEASE_DAY} ${ASLI_RELEASE_MONTH} ${ASLI_RELEASE_YEAR}")

# Info
set(DESCRIPTION "A Simple Lattice Infiller")
set(CMAKE_PROJECT_HOMEPAGE_URL "http://www.biomech.ulg.ac.be/ASLI")

# Build libraries as static (excluding dependencies related to parallelization)
set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so ${CMAKE_FIND_LIBRARY_SUFFIXES})
set(BUILD_SHARED_LIBS OFF)

# Graphical user interface default
option(ASLI_GUI "Compile the GUI" OFF)

# Documentation default
option(ASLI_DOC "Compile the documentation" OFF)

# Copy release dll's
option(ASLI_DLL "Copy required dll files" OFF)

# Mesher option defaults
option(CGAL_MESH "Include CGAL meshing option" ON)
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)
option(MMG_MESH "Include MMG meshing option" ON)
option(NATIVE "Add -march=native flag" OFF)

# Set some default paths
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules" ${CMAKE_ROOT}/modules)
if (NOT DEFINED EIGEN3_INCLUDE_DIR)
  set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libs/eigen)
endif()

# Compiler requirements
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.2.1)
    message(FATAL_ERROR "GCC version must be at least v10.2.1!")
  endif()
else()
  message(WARNING "You are using an unsupported compiler! Compilation has only been tested with GCC.")
endif()

# Keep Windows from breaking the min|max functions from the std library
if(WIN32) 
  add_definitions(-DNOMINMAX)
endif()


# Eigen
find_package(Eigen3 3.3.7 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# CGAL and its components
if(CGAL_MESH)
  add_definitions(-DCGAL_MESH)

  # CGAL
  find_package(CGAL 5.5.1 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR}/libs/CGAL) #Needs -frounding-math with gcc
  include_directories(${CGAL_DIR}/include)

  # Eigen (CGAL support)
  include(CGAL_Eigen_support)

  # Boost
  find_package(Boost 1.66 REQUIRED)

  # Compile for parallel use?
  if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3})
    add_definitions(-DCGAL_CONCURRENT_MESH_3)

    # TBB
    find_package(TBB REQUIRED)
    include(CGAL_TBB_support)
  endif()
endif()

# Mmg and its components
if(MMG_MESH)
  add_definitions(-DMMG_MESH)

  # MMG
  set(MMG_PREFIX mmg3d)
  set(MMG_PREFIX_DIR ${CMAKE_CURRENT_BINARY_DIR}/${MMG_PREFIX})
  set(MMG_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/${MMG_PREFIX})
  ExternalProject_Add(${MMG_PREFIX}
    PREFIX       ${MMG_PREFIX_DIR}
    SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/libs/mmg

    BUILD_ALWAYS OFF
    INSTALL_DIR  ${MMG_INSTALL_DIR}

    CMAKE_ARGS(-DCMAKE_BUILD_TYPE=Release -DBUILD=MMG3D -DLIBMMG3D_STATIC=ON
               -DLIBMMG3D_SHARED=OFF -DUSE_ELAS=OFF -DUSE_VTK=OFF -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)

    BUILD_COMMAND   make
    INSTALL_COMMAND make install
  )
  set(MMG3D_INCLUDE_DIRS "${MMG_INSTALL_DIR}/include")
  set(MMG3D_LIBRARIES "${MMG_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${MMG_PREFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")

  include_directories(${MMG3D_INCLUDE_DIRS})
  #add_library(MMG3D STATIC IMPORTED)
  #set_target_properties(MMG3D PROPERTIES IMPORTED_LOCATION "${MMG3D_LIBRARIES}")

  # MshMet
  set(MSHMET_PREFIX mshmet)
  set(MSHMET_PREFIX_DIR ${CMAKE_CURRENT_BINARY_DIR}/${MSHMET_PREFIX})
  set(MSHMET_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/${MSHMET_PREFIX})
  ExternalProject_Add(${MSHMET_PREFIX}
    PREFIX       ${MSHMET_PREFIX_DIR}
    SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/libs/AdaptTools

    BUILD_ALWAYS OFF
    BINARY_DIR ${MSHMET_INSTALL_DIR}
    INSTALL_DIR  ${MSHMET_INSTALL_DIR}

    CMAKE_ARGS(-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  )
  set(MSHMET_INCLUDE_DIRS "${MSHMET_INSTALL_DIR}/include")
  set(MSHMET_LIBRARIES "${MSHMET_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${MSHMET_PREFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")

  include_directories(${MSHMET_INCLUDE_DIRS})

  # Scotch
  find_package(SCOTCH)
  if(SCOTCH_FOUND)
    include_directories(PUBLIC ${SCOTCH_INCLUDE_DIRS})
  endif()

  # tetgen
  add_library(tet STATIC 
              ${CMAKE_CURRENT_SOURCE_DIR}/libs/tetgen/predicates.cxx 
              ${CMAKE_CURRENT_SOURCE_DIR}/libs/tetgen/tetgen.cxx
  )
  target_include_directories(tet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/tetgen)
  target_compile_definitions(tet PUBLIC TETLIBRARY) # -DTETLIBRARY: flag to compile tetgen as a library

endif()

# alglib
add_library(alg STATIC 
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/alglib/src/alglibinternal.cpp 
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/alglib/src/alglibmisc.cpp 
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/alglib/src/ap.cpp 
)
target_include_directories(alg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/alglib/src)

# mTT
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libs/mTT/include)

# yaml
file(GLOB yaml_SRC CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/libs/yaml/src/*.cpp) # Using file GLOB is not recomended!
add_library(yaml STATIC)
target_include_directories(yaml
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/yaml/include 
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/libs/yaml/src
)
target_sources(yaml PRIVATE ${yaml_SRC})
set_target_properties(yaml PROPERTIES CXX_STANDARD 11)

# Compile options for debuging
if(CMAKE_BUILD_TYPE MATCHES Debug)
  #add_compile_options(-W -Wall -Werror -Wextra -Wpedantic) # Better for W flags since it places them at top level
  #target_compile_options(-W -Wall -Werror -Wextra -Wpedantic)

  if(WIN32 OR MSYS) # Set to "large object" file format
  	add_compile_options(-Wa,-mbig-obj -O2)
  endif()
endif()

# Specify ASLI source files
set(ASLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/ASLI.cpp 
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/TrilinearInterpolation.cpp 
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Filter.cpp 
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Infill.cpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/icon.rc)

if(CGAL_MESH)
  set(ASLI_SOURCES ${ASLI_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/MeshCGAL.cpp) 
#                                   ${CMAKE_CURRENT_SOURCE_DIR}/src/OffsetCGAL.cpp) # TEMP!
endif()

if(MMG_MESH)
  set(ASLI_SOURCES ${ASLI_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/MeshMMG.cpp)
endif()

add_executable(ASLI ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${ASLI_SOURCES})

# Bring the headers into the build enviroment and set properties
target_include_directories(ASLI PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include) 
set_target_properties(ASLI PROPERTIES
                           RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin/$<0:>
                           CXX_STANDARD 17
)
if(MARCH_NATIVE)
  set_target_properties(ASLI PROPERTIES COMPILE_FLAGS "-march=native")
endif()

# Create entries for C++ files in "ASLI" routine
target_link_libraries(ASLI PUBLIC alg yaml)
if(NOT MSVC)
  target_link_libraries(ASLI PUBLIC stdc++fs)
endif()

if(CGAL_MESH)
  target_link_libraries(ASLI PUBLIC CGAL::CGAL CGAL::Eigen_support)
  if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
    target_link_libraries(ASLI PUBLIC TBB::tbb CGAL::TBB_support)
  endif()
endif()

if(MMG_MESH)
  add_dependencies(ASLI ${MMG_PREFIX} ${MSHMET_PREFIX})
  target_link_libraries(ASLI PUBLIC ${MMG3D_LIBRARIES} ${MSHMET_LIBRARIES} tet)
  if(SCOTCH_FOUND)
    target_link_libraries(ASLI PUBLIC ${SCOTCH_LIBRARIES} scotch)
  endif()
endif()

# Documentation
if(ASLI_DOC) # Automatic compiling of the documentation (unsupported feature)
  add_subdirectory(docs/manual docs/manual)

  add_custom_command(TARGET ASLI POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy
                       ${ASLI_MANUAL_BINARY_DIR}/manual.pdf
                       $<TARGET_FILE_DIR:ASLI>/docs/'ASLI [User Manual].pdf'
  )
else() # Copy pre-compiled documents
  add_custom_command(TARGET ASLI POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                       ${CMAKE_CURRENT_SOURCE_DIR}/docs/manual/manual.pdf
                       $<TARGET_FILE_DIR:ASLI>/docs/'ASLI [User Manual].pdf'
  )
endif()

# Copy config, license and demo files of ASLI to bin folder
add_custom_command(TARGET ASLI POST_BUILD
                          COMMAND ${CMAKE_COMMAND} -E copy_if_different
                            ${CMAKE_CURRENT_SOURCE_DIR}/config.yml
                            ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.*
                            ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
                            $<TARGET_FILE_DIR:ASLI>
                          COMMAND ${CMAKE_COMMAND} -E make_directory 
                            $<TARGET_FILE_DIR:ASLI>/inputs
                          COMMAND ${CMAKE_COMMAND} -E copy
                            ${CMAKE_CURRENT_SOURCE_DIR}/inputs/cube.*
                            $<TARGET_FILE_DIR:ASLI>/inputs
)

# The graphical user interface
if(ASLI_GUI)
  set(QASLI_PREFIX QASLI)
  set(QASLI_PREFIX_DIR ${CMAKE_CURRENT_BINARY_DIR}/${QASLI_PREFIX})
  set(QASLI_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

  ExternalProject_Add(${QASLI_PREFIX}
    PREFIX         ${QASLI_PREFIX_DIR}
    GIT_REPOSITORY https://github.com/tpms-lattice/QASLI.git
    GIT_TAG        main

    BUILD_ALWAYS OFF
    INSTALL_DIR ${QASLI_INSTALL_DIR}

    CMAKE_ARGS(-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
  )
endif()

# Copy release dll files (unsupported feature)
if(MSYS AND ASLI_DLL)
  if(CGAL_ACTIVATE_CONCURRENT_MESH_3)
    add_custom_command(TARGET ASLI POST_BUILD
  							COMMAND ${CMAKE_COMMAND} -E copy_if_different
  							${CMAKE_CURRENT_SOURCE_DIR}/dll/*.dll
  							$<TARGET_FILE_DIR:ASLI>
    )
  endif()

  if(ASLI_GUI)
    add_custom_command(TARGET ASLI POST_BUILD
   						      COMMAND ${CMAKE_COMMAND} -E copy_directory
   						      ${CMAKE_CURRENT_SOURCE_DIR}/dll/QASLI
   						      $<TARGET_FILE_DIR:ASLI>
    )
  endif()
endif()