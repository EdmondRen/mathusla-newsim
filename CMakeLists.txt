#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(Mu)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_BUILD_TYPE Debug)

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

#----------------------------------------------------------------------------
# Find CRY package
find_library(CRY NAMES CRY PATHS ${PROJECT_SOURCE_DIR}/cry_v1.7/lib NO_DEFAULT_PATH )
set(CRY_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/cry_v1.7/src)
set(CRY_LIBRARIES ${PROJECT_SOURCE_DIR}/cry_v1.7/lib/libCRY.a)
MESSAGE ("CRY installed at: " ${CRY_INCLUDE_DIRS} )



#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
                  ${CRY_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB_RECURSE  sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE  headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/include/*.hpp)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
## Main simulation executable
add_executable(simulation simulation.cc ${sources} ${headers})
target_link_libraries(simulation ${Geant4_LIBRARIES} ${CRY_LIBRARIES})

## Standalone CRY cosmic generator
add_executable(cry_ascii cry.cc ${headers})
target_link_libraries(cry_ascii ${CRY_LIBRARIES})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Mu. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

file(COPY macros DESTINATION ${PROJECT_BINARY_DIR})


# set(EXAMPLEMuA_SCRIPTS
#   exampleMu.out
#   exampleMu.in
#   gui.mac
#   init_vis.mac
#   plotHisto.C
#   plotNtuple.C
#   run1.mac
#   run2.mac
#   vis.mac
#   )

# foreach(_script ${EXAMPLEMuA_SCRIPTS})
#   configure_file(
#     ${PROJECT_SOURCE_DIR}/${_script}
#     ${PROJECT_BINARY_DIR}/${_script}
#     COPYONLY
#     )
# endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS simulation DESTINATION bin)
