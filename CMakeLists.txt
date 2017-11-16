cmake_minimum_required(VERSION 2.8)

# Define convenience macros

MACRO(HEADER_DIRECTORIES return_list start_dir)
  FILE(GLOB_RECURSE new_list "${start_dir}/*.h")
  SET(dir_list "")
  FOREACH(file_path ${new_list})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(dir_list ${dir_list} ${dir_path})
  ENDFOREACH()
  LIST(REMOVE_DUPLICATES dir_list)
  SET(${return_list} ${dir_list})
ENDMACRO()

# Project header

project("Electrostatic_Tools")

set(Electrostatic_Tools_VERSION_MAJOR 0)
set(Electrostatic_Tools_VERSION_MINOR 4)
set(Electrostatic_Tools_VERSION_PATCH 0)

set(Electrostatic_Tools_VERSION "${Electrostatic_Tools_VERSION_MAJOR}.${Electrostatic_Tools_VERSION_MINOR}.${Electrostatic_Tools_VERSION_PATCH}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Search for dependencies

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})
set(LIBS ${LIBS} ${OPENBABEL2_LIBRARIES})

find_package(Boost COMPONENTS program_options REQUIRED)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})
if(NOT EIGEN3_FOUND)
  message(WARNING "Eigen3 not found. mult_fitter and ep_fitter programs will not be built.")
  set(ET_NO_mult_fitter "YES")
  set(ET_NO_ep_fitter "YES")
endif(NOT EIGEN3_FOUND)

# Set includes 

configure_file("${CMAKE_SOURCE_DIR}/src/config.h.in" "${CMAKE_BINARY_DIR}/src/config.h")
include_directories("${CMAKE_BINARY_DIR}/src")

HEADER_DIRECTORIES(subdirs ${CMAKE_SOURCE_DIR}/src)
include_directories("${CMAKE_SOURCE_DIR}/src")
FOREACH(subdir ${subdirs})
    include_directories(${subdir})
ENDFOREACH()

# Set targets

file(GLOB sources "src/*.cpp")
foreach(file ${sources})
  get_filename_component (name ${file} NAME_WE)
  option(ET_NO_${name} "Do not build ${name} program" OFF)
  if(NOT ET_NO_${name})
    add_executable(${name} ${file})
    target_link_libraries(${name} ${LIBS})
    install(TARGETS ${name} DESTINATION bin)
  endif(NOT ET_NO_${name})
endforeach()

install(FILES "${CMAKE_SOURCE_DIR}/doc/manual.pdf" DESTINATION "share/doc/electrostatic-tools-${Electrostatic_Tools_VERSION}")
#install(DIRECTORY "${CMAKE_SOURCE_DIR}/doc/html" DESTINATION "share/doc/electrostatic-tools-${Electrostatic_Tools_VERSION}")
install(FILES "${CMAKE_SOURCE_DIR}/src/config/orient.rules" DESTINATION "share/electrostatic-tools/${Electrostatic_Tools_VERSION}")
install(FILES "${CMAKE_SOURCE_DIR}/src/config/placement.rules" DESTINATION "share/electrostatic-tools/${Electrostatic_Tools_VERSION}")

add_subdirectory(src/bindings)

# Source generation target

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Electrostatic Tools is a program package that aims to assist in developing of
next-generation molecular electrostatic models with an enhanced molecular electrostatic potential (MEP) anisotropy.")
SET(CPACK_PACKAGE_VENDOR "Eectrostatic Tools development team")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR ${Electrostatic_Tools_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${Electrostatic_Tools_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${Electrostatic_Tools_VERSION_PATCH})
SET(CPACK_SOURCE_IGNORE_FILES ".*toc;.*aux;.*dvi;.*4ct;.*4tc;.*idv;.*log;.*lg;.*xref;.*tmp")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "electrostatic-tools-${Electrostatic_Tools_VERSION}")

INCLUDE(CPack)


