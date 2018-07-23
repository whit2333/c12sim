cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(c12sim)

SET(C12SIM_MAJOR_VERSION 0)
SET(C12SIM_MINOR_VERSION 0)
SET(C12SIM_PATCH_VERSION 0)
set(C12SIM_VERSION ${C12SIM_MAJOR_VERSION}.${C12SIM_MINOR_VERSION}.${C12SIM_PATCH_VERSION})


SET(SHARE_DIR_NAME c12sim)
SET(MACRO_DIR      share/${SHARE_DIR_NAME})
SET(GEOMETRY_DIR   share/${SHARE_DIR_NAME}/geometry_files)

SET(C12SIM_MACRO_DIR ${CMAKE_INSTALL_PREFIX}/${MACRO_DIR} )
SET(C12SIM_GEOMETRY_DIR ${CMAKE_INSTALL_PREFIX}/${GEOMETRY_DIR} )

if(NOT C12SIM_DEFAULT_EG_DIR )
   SET(C12SIM_DEFAULT_EG_DIR ${CMAKE_INSTALL_PREFIX}/share/${SHARE_DIR_NAME}/data )
endif(NOT C12SIM_DEFAULT_EG_DIR )

if(NOT C12SIM_GDML_DIR )
   SET(C12SIM_GDML_DIR ${C12SIM_GEOMETRY_DIR})
endif(NOT C12SIM_GDML_DIR )

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

# CADMesh
find_package(cadmesh )
if(NOT cadmesh_FOUND)
  message("cadmesh NOT FOUND")
endif()

###################
# ClasDigi
find_package(ClasDigi REQUIRED)
include_directories(${ClasDigi_INCLUDE_DIRS})
set(C12SIM_LINK_LIBRARIES ${C12SIM_LINK_LIBRARIES} ${ClasDigi_LIBRARIES})
if(ClasDigi_DATA_DIR )
  SET(C12SIM_DATA_DIR ${ClasDigi_DATA_DIR})
endif(ClasDigi_DATA_DIR)
if(NOT C12SIM_DATA_DIR )
   SET(C12SIM_DATA_DIR ${CMAKE_INSTALL_PREFIX}/${SHARE_DIR_NAME}/noarch/data)
endif(NOT C12SIM_DATA_DIR )

###################
# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#set(ROOT_LIBS EG MathCore RIO Hist MathMore Graf Graf3d Gpad GX11 GeomPainter X3d Gviz3d RGL Eve EG Eve Foam GenVector Geom GLEW Gui HistPainter MathCore Matrix MLP Net Physics Spectrum Thread Tree TreeViewer EG)
set(ROOT_LIBS EG MathCore RIO Hist MathMore Graf Graf3d Gpad GeomPainter Gviz3d RGL Eve EG Eve Foam GenVector Geom Gui HistPainter MathCore Matrix MLP Net Physics Spectrum Thread Tree TreeViewer EG)
find_package(ROOT REQUIRED COMPONENTS ${ROOT_LIBS})
message(" ROOT_USE_FILE = ${ROOT_USE_FILE} " )
#include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})
set(C12SIM_LINK_LIBRARIES ${C12SIM_LINK_LIBRARIES} ${ROOT_LIBS})
message(" ROOT_LIBRARIES = ${ROOT_LIBRARIES} " )
message(" ROOT_LIBRARY_DIR = ${ROOT_LIBRARY_DIR} " )

#######################
# CADMesh (optional)
find_package(cadmesh)

link_directories(
   ${ClasDigi_LIBRARY_DIR}
   ${CLHEP_LIBRARY_DIR}
   ${ROOT_LIBRARY_DIR}
)
#
include_directories(
   ${PROJECT_SOURCE_DIR}/include
   ${PROJECT_BINARY_DIR})

# End of dependencies
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_BINARY_DIR}/C12SIM.h)

set(sources 
  #src/c12sim_help.cc
  src/helper_functions.cc
  src/B1ActionInitialization.cc
  src/B1Analysis.cc
  src/B1DetectorConstruction.cc
  src/B1DetectorMessenger.cc
  src/B1EventAction.cc
  src/B1OpticalPhysics.cc
  src/B1ParallelWorldConstruction.cc
  src/B1PrimaryGeneratorAction.cc
  src/B1RunAction.cc
  src/B1Run.cc
  src/B1SteppingAction.cc
  src/BeamlineDetectorGeometry.cc
  src/BeamOnTargetPrimaryGeneratorAction.cc
  src/BeamTestHit.cc
  src/BeamTestSD.cc
  src/C12MagneticField.cc
  src/DriftChamberDetectorGeometry.cc
  src/DriftChamberSensitiveDetector.cc
  src/ECDetectorGeometry.cc
  src/ECSensitiveDetector.cc
  src/FakeSD.cc
  src/FakeSDHit.cc
  src/GPSConfig.cc
  src/HTCCDetectorGeometry.cc
  src/HTCCSensitiveDetector.cc
  src/MicromegasVertexTrackerDetectorGeometry.cc
  src/PMTArraySensitiveDetector.cc
  src/RecoilChamberDetectorGeometry.cc
  src/RecoilChamberSensitiveDetector.cc
  src/RecoilHodoDetectorGeometry3.cc
  src/RecoilScintSensitiveDetector.cc
  src/ScintHit.cc
  src/SensitiveRegionDetector.cc
  src/SiliconVertexTrackerDetectorGeometry.cc
  src/SimplePrimaryGeneratorAction.cc
  src/SimulationManager.cc
  src/SimulationMessenger.cc
  src/SiPMHit.cc
  src/SiPMSD.cc
  src/SiTrackerSensitiveDetector.cc
  src/SiVertexTrackerDetectorGeometry.cc
  src/SolenoidDetectorGeometry.cc
  src/SteppingAction.cc
  src/TorusDetectorGeometry.cc
  src/TrackingAction.cc
  )
set(headers
  include/B1ActionInitialization.hh
  include/B1Analysis.hh
  include/B1DetectorConstruction.hh
  include/B1DetectorMessenger.hh
  include/B1EventAction.hh
  include/B1OpticalPhysics.hh
  include/B1ParallelWorldConstruction.hh
  include/B1PrimaryGeneratorAction.hh
  include/B1RunAction.hh
  include/B1Run.hh
  include/B1SteppingAction.hh
  include/BeamlineDetectorGeometry.h
  include/BeamOnTargetPrimaryGeneratorAction.hh
  include/BeamTestHit.hh
  include/BeamTestSD.hh
  include/C12MagneticField.h
  #include/c12sim_help.h
  #include/clipp.h
  include/dollar.hpp
  include/DriftChamberDetectorGeometry.h
  include/DriftChamberSensitiveDetector.h
  include/ECDetectorGeometry.h
  include/ECSensitiveDetector.h
  include/FakeSD.hh
  include/FakeSDHit.hh
  include/GPSConfig.hh
  include/HTCCDetectorGeometry.h
  include/HTCCSensitiveDetector.h
  include/MicromegasVertexTrackerDetectorGeometry.h
  include/PMTArraySensitiveDetector.h
  include/rang.hpp
  include/RecoilChamberDetectorGeometry.h
  include/RecoilChamberSensitiveDetector.h
  include/RecoilHodoDetectorGeometry3.h
  include/RecoilScintHit.hh
  include/RecoilScintSensitiveDetector.hh
  include/ScintHit.hh
  include/ScintTileDetectorGeometry.hh
  include/SensitiveRegionDetector.h
  include/SiliconVertexTrackerDetectorGeometry.h
  include/SimplePrimaryGeneratorAction.hh
  include/SimulationManager.h
  include/SimulationMessenger.h
  include/SiPMHit.hh
  include/SiPMSD.hh
  include/SiTrackerSensitiveDetector.h
  include/SiVertexTrackerDetectorGeometry.h
  include/SolenoidDetectorGeometry.h
  include/SteppingAction.hh
  include/TorusDetectorGeometry.h
  include/TrackingAction.hh
  )

#${PROJECT_BINARY_DIR}/C12SIM.h
include_directories(${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

add_executable(c12sim c12sim.cc ${sources} )

target_link_libraries(c12sim 
  PRIVATE ${Geant4_LIBRARIES}
  PRIVATE ${C12SIM_LINK_LIBRARIES} ${ClasDigi_LIBRARIES} ${ROOT_LIBRARIES}
  PRIVATE ${cadmesh_LIBRARIES}
  )
target_include_directories(c12sim 
  PRIVATE ${CADMESH_INCLUDE_DIRS} 
  )
target_compile_definitions(c12sim
  PRIVATE $<$<BOOL:${${cadmesh_FOUND}}>:USE_CADMESH>
)
#add_executable(c12sim c12sim.cc ${sources} ${headers})
#target_link_libraries(c12sim ${C12SIM_LINK_LIBRARIES} ${ClasDigi_LIBRARIES} ${ROOT_LIBRARIES})

target_compile_features(c12sim 
  PRIVATE cxx_std_14)

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
#install(TARGETS c12sim_wt DESTINATION bin)

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
  simple_eg_macro
  )
foreach(_script ${BASH_SCRIPT_FILES})
   configure_file("bin/${_script}.in" "${PROJECT_BINARY_DIR}/bin/${_script}" @ONLY)
   install(
      FILES         ${PROJECT_BINARY_DIR}/bin/${_script}
      PERMISSIONS  WORLD_READ WORLD_EXECUTE GROUP_EXECUTE  GROUP_READ OWNER_EXECUTE OWNER_WRITE OWNER_READ 
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

message("C12SIM_LINK_LIBRARIES = ${C12SIM_LINK_LIBRARIES} " )
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

