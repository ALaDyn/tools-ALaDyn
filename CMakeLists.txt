cmake_minimum_required (VERSION 3.1)
project(tools-ALaDyn)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_USE_MATH_DEFINES)
endif()

find_package(Boost COMPONENTS system filesystem REQUIRED)

set(JSONCONS_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/ext/jsoncons/src)

include_directories(${JSONCONS_INCLUDE_PATH} ${Boost_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src/)

add_executable(binary_analyzer                ${CMAKE_SOURCE_DIR}/binary_analyzer/binary_analyzer.cpp)
add_executable(diagspec_reader                ${CMAKE_SOURCE_DIR}/diagspec_reader/diagspec_reader.cpp)
add_executable(exponential_fit                ${CMAKE_SOURCE_DIR}/exponential_fit/exponential_fit.cpp)
add_executable(extract_nptot                  ${CMAKE_SOURCE_DIR}/extract_nptot/extract_nptot.cpp)
add_executable(fix_nptot                      ${CMAKE_SOURCE_DIR}/fix_nptot/fix_nptot.cpp)
add_executable(InGenUO                        ${CMAKE_SOURCE_DIR}/InGenUO/InGenUO.cpp)
target_link_libraries(InGenUO                 ${Boost_LIBRARIES}) 
add_executable(interpolate_scan_results       ${CMAKE_SOURCE_DIR}/interpolate_scan_results/interpolate_scan_results.cpp)
add_executable(lightweight_coredump_analyzer  ${CMAKE_SOURCE_DIR}/lightweight_coredump_analyzer/lightweight_coredump_analyzer.cpp)
add_executable(logaritmic_fit                 ${CMAKE_SOURCE_DIR}/logaritmic_fit/logaritmic_fit.cpp)
add_executable(merge_vtk                      ${CMAKE_SOURCE_DIR}/merge_vtk/merge_vtk.cpp)
add_executable(partdist_reader                ${CMAKE_SOURCE_DIR}/partdist_reader/partdist_reader.cpp)

add_executable(binary_decoder
  ${CMAKE_SOURCE_DIR}/binary_decoder/binary_decoder.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/binary_decoder.h
  ${CMAKE_SOURCE_DIR}/binary_decoder/binning.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/binning.h
  ${CMAKE_SOURCE_DIR}/binary_decoder/filter.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/filter.h
  ${CMAKE_SOURCE_DIR}/binary_decoder/grid_file_decoder.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/grid_file_decoder.h
  ${CMAKE_SOURCE_DIR}/binary_decoder/parameters.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/parameters.h
  ${CMAKE_SOURCE_DIR}/binary_decoder/phasespace_file_decoder.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/phasespace_file_decoder.h
  ${CMAKE_SOURCE_DIR}/binary_decoder/swap_tools.cpp
  ${CMAKE_SOURCE_DIR}/binary_decoder/swap_tools.h
)

install (TARGETS 
  binary_analyzer
  binary_decoder
  diagspec_reader
  exponential_fit
  extract_nptot
  fix_nptot
  InGenUO
  interpolate_scan_results
  lightweight_coredump_analyzer
  logaritmic_fit
  merge_vtk
  partdist_reader
  DESTINATION ${CMAKE_SOURCE_DIR}/bin/
)

