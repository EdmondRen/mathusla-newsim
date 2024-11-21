# mathusla-newsim


# Standalone generator: CRY

```
CRY cosmic generator with text output
Usage:
  CRY cosmic generator \[OPTION...]

  -h, --help          Print help   
  -s, --setup arg     setup file name (default:   
                      ../macros/generators/cry_all.file)   
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

    ./cry_ascii -c 2112 -e 2600 -n 1000000