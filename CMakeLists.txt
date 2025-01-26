#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(Mu)

## Compiler settings
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)


## Build type: Debug or Release
# set (CMAKE_BUILD_TYPE Release)
set (CMAKE_BUILD_TYPE Debug)
option(PROFILE_WITH_GPROF "Enable profiling with gprof" OFF)

    
## Enable profiling for Debug build
# Run the simulation, then use `gprof my_program gmon.out > analysis.txt` to check the result
# if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Profile")
if(PROFILE_WITH_GPROF)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    add_compile_options(-pg)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

#----------------------------------------------------------------------------
# Find packages

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Find CRY package
find_library(CRY NAMES CRY PATHS ${PROJECT_SOURCE_DIR}/external/cry_v1.7/lib NO_DEFAULT_PATH )
set(CRY_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external/cry_v1.7/src)
set(CRY_LIBRARIES ${PROJECT_SOURCE_DIR}/external/cry_v1.7/lib/libCRY.a)
MESSAGE ("CRY installed at: " ${CRY_INCLUDE_DIRS} )

# ROOT package. Only used for digitizer and metadata
find_package(ROOT REQUIRED)

# Eigen3
find_package (Eigen3 REQUIRED NO_MODULE)


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE} )
include(${ROOT_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${CRY_INCLUDE_DIRS}
                    ${EIGEN3_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB_RECURSE  sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE  headers ${PROJECT_SOURCE_DIR}/include/*.hh
                           ${PROJECT_SOURCE_DIR}/include/*.hpp)

file(GLOB  sources_tracker ${PROJECT_SOURCE_DIR}/tracker/*.cc)
file(GLOB  headers_tracker ${PROJECT_SOURCE_DIR}/tracker/*.hh)                           

#----------------------------------------------------------------------------
# Make library and executables

# Make a shared library for simulation
add_library(mu_sim SHARED ${sources})
target_include_directories(mu_sim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(mu_sim PUBLIC
    ${Geant4_LIBRARIES} 
    ${ROOT_LIBRARIES})

# Make a shared library for tracker
add_library(mu_recon SHARED ${sources_tracker})
target_include_directories(mu_recon PUBLIC
    ${PROJECT_SOURCE_DIR}/tracker/
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_SOURCE_DIR}/include>)
target_link_libraries(mu_recon PUBLIC
    ${ROOT_LIBRARIES})    

## Simulation main entrance
add_executable(simulation simulation.cc)
target_link_libraries(simulation PUBLIC mu_sim ${CRY_LIBRARIES})


## Digitizer
add_executable(digitizer digitizer.cc)
target_link_libraries(digitizer PUBLIC mu_sim ${CRY_LIBRARIES})

## Standalone CRY cosmic generator
add_executable(cry cry.cc ${headers})
target_link_libraries(cry PUBLIC mu_sim ${CRY_LIBRARIES})

## Standalone parma cosmic generator
add_executable(parma parma.cc)
target_link_libraries(parma PUBLIC mu_sim ${CRY_LIBRARIES})

## Tracker tests
add_executable(tracker_test1 ${PROJECT_SOURCE_DIR}/tracker/tests/test_kalman_model.cc)
target_link_libraries(tracker_test1 mu_recon)
add_executable(tracker_test2 ${PROJECT_SOURCE_DIR}/tracker/tests/test_track_finder.cc)
target_link_libraries(tracker_test2 mu_recon)


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Mu. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

# file(COPY macros DESTINATION ${PROJECT_BINARY_DIR})

# set(EXAMPLEMuA_SCRIPTS
#   exampleMu.out
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
# install(TARGETS cry DESTINATION bin)
