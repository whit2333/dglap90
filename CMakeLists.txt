# ==============================================================================
# CMakeLists.txt
#
# by Adam Freese <maxwellsdemon137@gmail.com> copyright 2017
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This file is part of dglap90.
#
# dglap90 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dglap90 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with dglap90.  If not, see <http://www.gnu.org/licenses/>.

# ==============================================================================
# //////////////////////////////////////////////////////////////////////////////
# Configuration
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ==============================================================================

cmake_minimum_required (VERSION 2.6)
find_program(CMAKE_Fortran_COMPILER $ENV{FC} ifort gfortran PATH)

project(dglap)

enable_language(Fortran)

# ==============================================================================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
# ==============================================================================
get_filename_component(FCN ${CMAKE_Fortran_COMPILER} NAME)
if(FCN MATCHES "gfortran.*")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fopenmp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -fbounds-check -fcheck=all")
elseif(FCN MATCHES "ifort.*")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fopenmp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback")
else()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g")
endif()
# ==============================================================================
if(CMAKE_BUILD_TYPE MATCHES "RELEASE")
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES "DEBUG")
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
endif()
# ==============================================================================
site_name(computer_name)
message("cmake ver : " ${CMAKE_VERSION}                  )
message("project   : " ${PROJECT_NAME}                   )
message("hostname  : " ${computer_name}                  )
message("system    : " ${CMAKE_SYSTEM}                   )
message("arch      : " ${CMAKE_HOST_SYSTEM_PROCESSOR}    )
message("compiler  : " ${CMAKE_Fortran_COMPILER}         )
message("command   : " ${FCN} " " ${CMAKE_Fortran_FLAGS} )
message("build     : " ${CMAKE_BUILD_TYPE}               )
message("sources   : " ${CMAKE_CURRENT_SOURCE_DIR}       )
# ==============================================================================
set(LIBRARY_OUTPUT_PATH            "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH         "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
# ==============================================================================

# ==============================================================================
# //////////////////////////////////////////////////////////////////////////////
# Libraries
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ==============================================================================

set(SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

set(LIBS "basics" "dataframes" "qcd" "dglap")
foreach(LIB ${LIBS})
  add_library(${LIB} SHARED ${SOURCE_PATH}/${LIB}.f90)
endforeach()

# Everything else depends on basics
foreach(LIB ${LIBS})
  if(NOT ${LIB} MATCHES "basics")
    target_link_libraries(${LIB} PUBLIC basics)
  endif()
endforeach()

# dglap depends on everything else
foreach(LIB ${LIBS})
  if(NOT ${LIB} MATCHES "dglap")
    target_link_libraries(dglap PUBLIC ${LIB})
  endif()
endforeach()

# ==============================================================================
# //////////////////////////////////////////////////////////////////////////////
# Test executable
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ==============================================================================

add_executable(test.bin ${SOURCE_PATH}/tester.f90)

foreach(LIB ${LIBS})
  target_link_libraries(test.bin ${LIB})
endforeach()

# ==============================================================================
# //////////////////////////////////////////////////////////////////////////////
# Copy PDF grids and plotter from source to build directory
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ==============================================================================

set( PDFGRIDS "cj15_Q2_5_GeV2.csv" "cj15_Q2_1000_GeV2.csv" )
foreach( GRID ${PDFGRIDS} )
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pdf_grids/${GRID}
    ${CMAKE_BINARY_DIR}/pdf_grids/${GRID} COPYONLY)
  message( "Resource file : ${GRID}" )
endforeach()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plots.py.in
  ${CMAKE_BINARY_DIR}/plots.py)
