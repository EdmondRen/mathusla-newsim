#ifndef Uoft1_Builder_h
#define Uoft1_Builder_h 1

// STD libs
#include <any>


// G4 include
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// User include
#include "geometry/_GeoBuilder.hh"


/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

namespace MuGeoBuilder
{


    // Geometry Builder Class 
    // This class takes care of detector construction 
    // as well as assigning sensitive detector to correct volumes
    class Uoft1_Builder : public Builder
    {
    public:
        Uoft1_Builder(const std::string &detector_name);

        // Core function 1:
        // Construct physics volume, should return world PV
        G4VPhysicalVolume *Construct();

        // Core function 2:
        // Set the sensitive detector for this geometry

        // Helper functions:
        void DefineMaterials();
        void DefineGeometry();
        G4VPhysicalVolume* ConstructLayer();
        G4VPhysicalVolume* ConstructModule();
        G4VPhysicalVolume* ConstructDetector();
        G4VPhysicalVolume* ConstructEnvironment();

        // std::unordered_map<std::string, std::any> GeometryConfig;
    };

} //namespace MuGeoBuilder

#endif

