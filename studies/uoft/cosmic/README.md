# Cosmic simulation for teststand


`cry_all_cms_0m.conf` file contains the CRY cosmic generator settings. subboxLength is the length of a square area that will be sampled. 

`g4config_cry_all_mathusla40_abstime.mac` contains Geant4 simulation settings. 
"/gen/cry/box 2. 2. 2. m" creates a box of 2x2x2 m. This is used to preselect the particles.
"/gen/cry/offset 0. 0. 2. m" offsets the box by 2m upwards so that it is above ground. 


Command example:

(in ./build folder)

```bash
./simulation -r 1 -s 1 \
    -m ../studies/uoft/cosmic/g4config_cry_all_mathusla40_abstime.mac,events,1000 \
    -o $YOUR_DATA_FOLDER
```