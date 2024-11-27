# mathusla-newsim

## General controls

### About random number generator

This simulation is built around one single static generator, G4Random::RanecuEngine

Notes:
* Do not use other engines. This is the only engine that can be restored to certain status with just 3 integers.
* If you want to use the engine in other parts of the code, use CLHEP::HepRandom::getTheEngine()
* The engien status after the generator is saved in the result as \[Seed_init, Seed_0, Seed_1\]
  * Seed are saved as integer in the ROOT file, but they need to be converted to unsigned int before using it.
  * Example: `*reinterpret_cast<unsigned int*>(&seed)`

## Event generator

### For user

### For developer

## Geometry

### For user

### For developer


## Reproducing an event based on saved status



# Digitizer

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