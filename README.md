# MATHUSLA simulation V2

## 1. General controls

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

### 2.4. FileReader generator (filereader)

### 2.5. Generator for recreating previous events (recreate)

This generator is able to re-generate previous events, based on the saved random number generator status (two integers) and primary vertex information in the ROOT file.
It takes in a list of event records. Each line contains a {filename, entry} pair, where filename is the ROOT file to take the information, and entry is the event entry in the ROOT file.



## 3. Geometry








# Event Digitizer

Digitizer takes the simulation output ROOT file, turns the simulation truth into measurable quantities (x, y, z, t) along with some other information for debugging. 

## Digitizer output format

There are two ROOT tuples in the digitizer output ROOT file: "digi" and "metadata".

"metadata" tree contains general information, while "digi" contains event-wise information.

| Key             | dtype | Comment                                                                                                                                                                                                                                                                                                                                                                                 |
|-----------------|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Digi_{x,y,z,t}  | float | Coordinate and time of digis                                                                                                                                                                                                                                                                                                                                                            |
| Digi_edep       | float | Energy deposition                                                                                                                                                                                                                                                                                                                                                                       |
| Digi_trackID    | int32 | Track ID, index of the track in Geant4                                                                                                                                                                                                                                                                                                                                                  |
| Digi_pdgID      | int32 | Particle PDG ID                                                                                                                                                                                                                                                                                                                                                                         |
| Digi_detectorID | int64 | Detector ID, in the format of AAABBBCCCDDDDD. A, B, C, Ds are the copy number of detector, tower-module, layer, bar.                                                                                                                                                                                                                                                                    |
| Digi_type       | int32 | Indicates which generator the digi is from. 0: gun, 1: cry, 2: parma, -1: noise                                                                                                                                                                                                                                                                                                         |
| Digi_hitInds    | int32 | Indices of the truth hits of each digi. Separated by -1.                                                                                                                                                                                                                                                                                                                                |
| Digi_direction  | int32 | Direction of the bar coded into the last three digits of this number. Each digit indicates the current direction of the bar. The hundreds place is for the original x direction, the tens place for the original y direction, and the ones place for x. For example, 201 means the x direction (along the bar) is now pointing to z, y direction is now pointing to x, and z is now y.  |

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
