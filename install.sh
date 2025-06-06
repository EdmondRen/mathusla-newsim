#!/bin/bash

# MATHUSLA MU Detector Simulation : Build Script
#
# Copyright 2018 Brandon Gomes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

CONFIG_CMAKE=0
CLEAN_BUILD=0
RUN_BUILD=0
OTHER_ARGS=""
EIGEN3_USER="$HOME/geant_projects/lib/share/eigen3/cmake"

while [[ $# -gt 0 ]]; do
    arg="$1"
    case $arg in
        --cmake | cmake) CONFIG_CMAKE=1
            shift;;
        --clean | clean) CLEAN_BUILD=1
            shift;;
        --run | run) RUN_BUILD=1
            shift;;
        --help | help | -h)
            echo -e "usage: ./install [--clean] [--cmake] [--run [args]]"
            echo -e "  clean      : rebuild all source files"
            echo -e "  cmake      : rebuild with CMake"
            echo -e "  run [args] : run executable after building with \"args\""
            echo -e "  help       : open this help screen"
            exit;;
        *) OTHER_ARGS="$OTHER_ARGS $arg"
            shift;;
    esac
done




# Build CRY first
echo "------------------------------------------------"
echo "  Building CRY"
cd external/cry_v1.7
make clean
make -j8
cd ../..
echo "  Finished building CRY"
echo "------------------------------------------------\n"


if [[ "$CONFIG_CMAKE" -eq 1 && "$CLEAN_BUILD" -eq 1 ]]; then
    rm -rf build;
    rm simulation;
fi

mkdir -p build
cd build

if [[ "$CONFIG_CMAKE" -eq 1 ]]; then cmake -DCMAKE_INSTALL_PREFIX=`realpath ..` .. -DEigen3_DIR=$EIGEN3_USER;   fi;
if [[ "$CLEAN_BUILD"  -eq 1 ]]; then make clean; fi;
pwd

echo "------------------------------------------------"
echo "  Building MATHUSLA simulations"
make -j8
# make install -j8
echo "  Finished building MATHUSLA simulations"
echo "------------------------------------------------"

if [[ ! "$?" -eq 0 ]]; then
    echo -e "\nBuild Failed!\n"
    exit 1
else
    # mv simulation ../simulation
    if [[ "$RUN_BUILD" -eq 1 ]]; then cd ..; ./simulation $OTHER_ARGS; fi;
    exit $?
fi
                                                                                               