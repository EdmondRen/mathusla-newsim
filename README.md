# MATHUSLA simulation V2

This repository contains the simulation and reconstruction code for MATHUSLA experiment. 

## Install

### Dependencies

There are three external libraries needed to complie this simulation code

1. GEANT4 10.07
2. ROOT 6.28/10
3. Eigen 3: https://eigen.tuxfamily.org/index.php?title=Main_Page

GEANT4 and ROOT have detailed instruction on the website using CMake build system. They should be compliled with C++ 17 standard as it is used by this package. 

Eigen is a header-only library. It doesn't matter where you put the downloaded file, but you should still run cmake to build the automate detection files for cmake. After that a Eigen3Config.cmake file is generate, which is usually at lib/share/eigen3/cmake. This path needs to be saved in install.sh as $EIGEN3_USER

On cedar, you can load the requried packages with command:

    module load StdEnv/2020 gcc/9.3.0 qt/5.12.8 root/6.26.06  eigen/3.3.7 geant4/10.7.3 geant4-data/10.7.3

### How to build

There is a shell script `install.sh` to automate the build process. For the first time, run `./install --cmake`. Aftwards only `./install` is needed.

    usage: ./install [--clean] [--cmake] [--run [args]]
      clean      : rebuild all source files
      cmake      : rebuild starting with CMake configuration
      help       : open this help screen

### Binaries

There are six executables saved in `.\build\` folder after the build finishes.

| Name | Description|
|---|---|
| `simulation` | main simulation executable |
| `digitizer`  | convert simulation truth to measurable quantities in detector |
| `tracker`    | reconstruct track and vertex from the digitized events | 
| `param` | Standalone cosmic ray event generator|
| `cry` | Standalone cosmic ray event generator|

# `simulation`

## 1. General controls

Options for simulation:

    Usage:
      ./simulation [OPTION...]

      -h, --help           Print help
      -D, --debug          Enable debugging
      -m, --macro arg      Macro name, followed by possible parameters for 
                          macro seperated by commas. For example, 
                          -m=run1.mac,Ek,10,theta,20
      -g, --generator arg  Generator, one of <gun/range/parma/cry/filereader> 
                          (default: gun)
      -d, --detector arg   Detector, one of <math40/uoft1> (default: uoft1)
      -e, --export arg     Export directory for geometry and other information. 
                          Default is ./export/ (default: export)
      -o, --output arg     Output directory. Default is ./data/ (default: data)
      -r, --run arg        Run number (default: 0)
      -s, --seed arg       Seed of random number generator, a positive integer. 
                          Default to -1 for random seed. Events in the same 
                          run share the same seed. (default: -1)
      -S, --session arg    Session name (default: MathuslaSim)
      -a, --all_steps      [Warning] Use cautiously. Write all steps into ROOT 
                          file. This usually takes a lot of storage space and 
                          more time on IO.
      -v, --save_vrml      [Warning] Use cautiously. Write visualization of 
                          each event into VRML file. This usually takes a lot 
                          of storage space and more time on IO.
      -t, --threads arg    [NOT IMPLEMENTED] Number of threads. Only works with 
                          single (1) thread at this moment. (default: 1)

### About random number generator

This simulation is built around one single static generator, G4Random::RanecuEngine

Notes:
* Do not use other engines. This is the only engine that can be restored to certain status with just 3 integers.
* If you want to use the engine in other parts of the code, use CLHEP::HepRandom::getTheEngine()
* The engien status after the generator is saved in the result as \[Seed_init, Seed_0, Seed_1\]
  * Seed are saved as integer in the ROOT file, but they need to be converted to unsigned int before using it.
  * Example: `*reinterpret_cast<unsigned int*>(&seed)`

## 2. Event generator

### 2.1. Geant4 prticle gun (gun)

### 2.2. PARMA cosmic generator (parma)

### 2.3. CRY cosmic generator (cry)

### 2.5. Generator for recreating previous events (recreate)

This generator is able to re-generate previous events, based on the saved random number generator status (two integers) and primary vertex information in the ROOT file. It can also be used to creat events by simply providing the {position, momentum, pdg_number} of each particle. 

It takes in a list of event records. Each line contains a {filename, entry} pair, where filename is the ROOT file to take the information, and entry is the event entry in the ROOT file.



## 3. Geometry








# Event Digitizer

Digitizer takes the simulation output ROOT file, turns the simulation truth into measurable quantities (x, y, z, t) along with some other information for debugging. 

    Usage:
      CRY cosmic generator [OPTION...] positional parameters

      -h, --help                  Print help
      -s, --seed arg              Seed for random number generator (default: 
                                  -1)
      -t, --time_resolution arg   Coincidence time resolution [ns]. (default: 
                                  1)
      -T, --time_limit arg        Time limit [ns] (default: 20)
      -E, --energy_threshold arg  Energy threshold for a digi [MeV] (default: 
                                  0.65)
      -n, --noise arg             Noise rate [avg number per file]. Set to -1 
                                  to disable (default). (default: -1)
      -p, --print_progress arg    Print progress every `p` events (default: 1)
      -w, --window arg            Noise window [ns]

## Digitizer output format

There are two ROOT tuples in the digitizer output ROOT file: "digi" and "metadata".

"metadata" tree contains general information, while "digi" contains event-wise information.

| Key             | dtype | Comment                                                                                                                                                                                                                                                                                                                                                                                 |
|-----------------|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Digi_{x,y,z,t}  | float | Coordinate [mm] and time [ns] of digis                                                                                                                                                                                                                                                                                                                                                            |
| Digi_edep       | float | Energy deposition                                                                                                                                                                                                                                                                                                                                                                       |
| Digi_trackID    | int32 | Track ID, index of the track in Geant4                                                                                                                                                                                                                                                                                                                                                  |
| Digi_pdgID      | int32 | Particle PDG ID                                                                                                                                                                                                                                                                                                                                                                         |
| Digi_detectorID | int64 | Detector ID, in the format of AAABBBCCCDDDDD. A, B, C, Ds are the copy number of detector, tower-module, layer, bar.                                                                                                                                                                                                                                                                    |
| Digi_type       | int32 | Indicates which generator the digi is from. 0: gun, 1: cry, 2: parma, -1: noise                                                                                                                                                                                                                                                                                                         |
| Digi_hitInds    | int32 | Indices of the truth hits of each digi. Separated by -1.                                                                                                                                                                                                                                                                                                                                |
| Digi_direction  | int32 | Direction of the bar coded into the last three digits of this number. Each digit indicates the current direction of the bar. The hundreds place is for the original x direction, the tens place for the original y direction, and the ones place for x. For example, 201 means the x direction (along the bar) is now pointing to z, y direction is now pointing to x, and z is now y.  |


# Python helper functions

Python helper function are located in `/python` folder. The library can be installed with the following commands:

```bash
cd python
pip install -e . --user
```

# Standalone generator: CRY

```
CRY cosmic generator with text output
Usage:
  CRY cosmic generator \[OPTION...]

  -h, --help          Print help   
  -s, --setup arg     setup file name (default:   
                      ../macros/generators/cry_all.conf)   
  -n, --nevents arg   Number of events (default: 1000)   
  -e, --ekin arg      Kinetic energy cut. Positive value means above, and   
                      negative means below. Set to 0 to disable (default:   
                      0)  
  -c, --contains arg  Select only events that contains specific PDG-ids,   
                      separated by commas. Eg, --contains=13,-13. All   
                      supported: e(+-11), mu(+-13), neutron(2112),   
                      proton(2212), gamma(22), pion(+-211,111),   
                      kaon(+-321,311,310,130). Set to 0 to disable   
                      (default: 0)  
```                    

Example: Generate 1 million events, save the ones that contains neutrons above 2600 MeV.

    ./cry -c 2112 -e 2600 -n 1000



./simulation -m ../macros/run_cry.mac -s 1
./cry -c 2212 -n 10000 -s ../macros/generators/cry_proton.conf
./cry -n 10000 -s ../macros/generators/cry_proton.conf


# Standalone generator: PARMA


# Credits

## json config parser

https://github.com/taocpp/config/tree/main