cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(c12sim)

SET(C12SIM_MAJOR_VERSION 0)
SET(C12SIM_MINOR_VERSION 0)
SET(C12SIM_PATCH_VERSION 0)
set(C12SIM_VERSION ${C12SIM_MAJOR_VERSION}.${C12SIM_MINOR_VERSION}.${C12SIM_PATCH_VERSION})

SET(C12SIM_MACRO_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME} )

if(NOT C12SIM_DATA_DIR )
   SET(C12SIM_DATA_DIR $ENV{JLAB_ROOT}/noarch/data)
endif(NOT C12SIM_DATA_DIR )

if(NOT C12SIM_DEFAULT_EG_DIR )
   SET(C12SIM_DEFAULT_EG_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/data )
endif(NOT C12SIM_DEFAULT_EG_DIR )

if(NOT C12SIM_GDML_DIR )
   SET(C12SIM_GDML_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/data/gdml )
endif(NOT C12SIM_GDML_DIR )


#set(CMAKE_CXX_FLAGS "-g3 ${CMAKE_CXX_FLAGS} ")
SET(C12SIM_LINK_LIBRARIES)

#----------------------------------------------------------------------------
# Dependencies
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
set(C12SIM_LINK_LIBRARIES ${C12SIM_LINK_LIBRARIES} ${Geant4_LIBRARIES})

###################
# ClasDigi
find_package(ClasDigi REQUIRED)
#include(${ClasDigi_USE_FILE})
include_directories(${ClasDigi_INCLUDE_DIRS})
set(C12SIM_LINK_LIBRARIES ${C12SIM_LINK_LIBRARIES} ${ClasDigi_LIBRARIES})

###################
# CLHEP
set (CLHEP_CONFIG_DEBUG 1)
find_package(CLHEP     REQUIRED)
#  CLHEP_INCLUDE_DIRS - include directories for CLHEP
#  CLHEP_DEFINITIONS  - compile definitions needed to use CLHEP
#  CLHEP_LIBRARIES    - libraries to link against
include_directories(${CLHEP_INCLUDE_DIRS})
add_definitions(${CLHEP_DEFINITIONS})
#message(" CLHEP_DEFINITIONS = ${CLHEP_DEFINITIONS} " )
set(C12SIM_LINK_LIBRARIES ${C12SIM_LINK_LIBRARIES} ${CLHEP_LIBRARIES})

###################
# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(ROOT_LIBS EG MathCore RIO Hist MathMore Graf Graf3d Gpad GX11 GeomPainter X3d Gviz3d RGL Eve EG Eve Foam GenVector Geom GLEW Gui HistPainter MathCore Matrix MLP Net Physics Spectrum Thread Tree TreeViewer EG)
find_package(ROOT REQUIRED COMPONENTS ${ROOT_LIBS})
message(" ROOT_USE_FILE = ${ROOT_USE_FILE} " )
#include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})
set(C12SIM_LINK_LIBRARIES ${C12SIM_LINK_LIBRARIES} ${ROOT_LIBS})
message(" ROOT_LIBRARIES = ${ROOT_LIBRARIES} " )
message(" ROOT_LIBRARY_DIR = ${ROOT_LIBRARY_DIR} " )


link_directories(
   ${ClasDigi_LIBRARY_DIR}
   ${CLHEP_LIBRARY_DIR}
   ${ROOT_LIBRARY_DIR}
)

include_directories(
   ${PROJECT_SOURCE_DIR}/include
   ${PROJECT_BINARY_DIR})

###################
message(" C12SIM_LINK_LIBRARIES = ${C12SIM_LINK_LIBRARIES} " )
# End of dependencies
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_BINARY_DIR}/C12SIM.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(c12sim c12sim.cc ${sources} ${headers})
target_link_libraries(c12sim ${C12SIM_LINK_LIBRARIES} ${ClasDigi_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(C12SIM DEPENDS c12sim)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MACROB1_SCRIPTS
  mac/init_vis.mac
  mac/run1.mac
  mac/run2.mac
  mac/vis.mac
  mac/vis_parallel.mac
  mac/fastMC.mac
  )
foreach(_script ${MACROB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
install(FILES ${MACROB1_SCRIPTS} DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}" )


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS c12sim DESTINATION bin)

# ----------------------------------------------------------------------------
# Configured files 
# These files end with ".in" which is removed when installed
set(EXAMPLEB1_CONF_FILES
  c12sim.h
  init_default.mac
  #vis.mac
  )
foreach(_script ${EXAMPLEB1_CONF_FILES})
  configure_file("${_script}.in" "${PROJECT_BINARY_DIR}/${_script}" @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/${_script} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME} )
endforeach()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/data DESTINATION ${C12SIM_DEFAULT_EG_DIR}/.. )

#------------------
# bash scripts
set(BASH_SCRIPT_FILES
  fast-c12sim
  )
foreach(_script ${BASH_SCRIPT_FILES})
   configure_file("bin/${_script}.in" "${PROJECT_BINARY_DIR}/bin/${_script}" @ONLY)
   install(
      FILES         ${PROJECT_BINARY_DIR}/bin/${_script}
      PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ 
      DESTINATION bin )
endforeach()

#-------------------

## Add all targets to the build-tree export set
export(TARGETS c12sim FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
#
## Export the package for use from the build-tree
## (this registers the build-tree with a global CMake-registry)
#export(PACKAGE ${PROJECT_NAME})
#
#set(CONF_LIBRARIES ${PROJECT_NAME})
#
#set(CONF_EXECUTABLES insane-config)
#
#set(CONF_INCLUDE_DIRS 
#   ${CMAKE_INSTALL_PREFIX}/include
#   ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
#   )
#
#set(CONF_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib)
#
##configure_file(${PROJECT_NAME}Config.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake" @ONLY)
#
#configure_file(${PROJECT_NAME}ConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake" @ONLY)
#
## Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
#install(FILES
#  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake"
#  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
#  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
#
## Install the export set for use with the install-tree
#install(EXPORT ${PROJECT_NAME}Targets DESTINATION
#  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

message("C12SIM_MAJOR_VERSION    ${C12SIM_MAJOR_VERSION}")
message("C12SIM_MINOR_VERSION    ${C12SIM_MINOR_VERSION}")
message("C12SIM_PATCH_VERSION    ${C12SIM_PATCH_VERSION}")
message("C12SIM_MACRO_DIR        ${C12SIM_MACRO_DIR}")
message("C12SIM_DATA_DIR         ${C12SIM_DATA_DIR}")
message("C12SIM_DEFAULT_EG_DIR   ${C12SIM_DEFAULT_EG_DIR}")
message("C12SIM_CXX_CFLAGS       ${C12SIM_CXX_CFLAGS}")
message("C12SIM_CXX_LDFLAGS      ${C12SIM_CXX_LDFLAGS}")
message("C12SIM_CXX_LIBS         ${C12SIM_CXX_LIBS}")
message("C12SIM_CXX_INC_DIR      ${C12SIM_CXX_INC_DIR}")
message("C12SIM_PREFIX           ${CMAKE_INSTALL_PREFIX}")

