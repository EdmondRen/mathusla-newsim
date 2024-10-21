#ifndef GeoBuilder_HH
#define GeoBuilder_HH 1

#include "G4Material.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GenericMessenger.hh"
#include <G4UIcommand.hh>


namespace MuGeoBuilder
{


    // Geometry Builder Class 
    // This class takes care of detector construction 
    // as well as assigning sensitive detector to correct volumes
    class Builder
    {
    public:
        Builder();

        // Core function 1:
        // Construct physics volume, should return world PV
        virtual G4VPhysicalVolume *Construct();
        G4VPhysicalVolume *worldPV;

        // Core function 2:
        // Set the sensitive detector for this geometry


        // Messenger related
        G4GenericMessenger *fMessenger;
        std::string DetectorName;

    private:
    };


    namespace Material
    {
        extern G4Element *H;
        extern G4Element *C;
        extern G4Element *N;
        extern G4Element *O;
        extern G4Element *F;
        extern G4Element *Al;
        extern G4Element *Si;
        extern G4Element *S;
        extern G4Element *Ar;
        extern G4Element *Ca;
        extern G4Material *Air;
        extern G4Material *Aluminum;
        extern G4Material *Bakelite;
        extern G4Material *Copper;
        extern G4Material *Concrete;
        extern G4Material *Iron;
        extern G4Material *PolystyreneFoam;
        extern G4Material *Polyvinyltoluene;
        extern G4Material *PlasticScintillator;
        extern G4Material *LiquidArgon;
        extern G4Material *Vacuum;
        
    } /* namespace Material */    

} // namespace MuGeoBuilder


#endif
