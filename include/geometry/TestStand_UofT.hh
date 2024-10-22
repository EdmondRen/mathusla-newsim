#ifndef Uoft1_Builder_h
#define Uoft1_Builder_h 1

// STD libs
#include <any>

// G4 include
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// User include
#include "geometry/_GeoBuilder.hh"

namespace MuGeoBuilder
{
  // Geometry Builder Class
  // This class takes care of detector construction
  // as well as assigning sensitive detector to correct volumes
  class Uoft1_Builder : public Builder
  {
  public:
    Uoft1_Builder();

    // Core function 1:
    // Construct physics volume, should return world PV
    G4VPhysicalVolume *Construct() override;
    std::vector<G4VPhysicalVolume *> allSensitiveDetectors;
    G4LogicalVolume *worldLV;
    G4LogicalVolume *detectorLV;
    G4LogicalVolume *environmentLV;

    // Core function 2:
    // Set the sensitive detector for this geometry

    // Helper functions:
    // void DefineMaterials();
    // void DefineGeometry();
    G4LogicalVolume *ConstructLayer(G4LogicalVolume *parentVolume);
    G4LogicalVolume *ConstructModule(G4LogicalVolume *parentVolume);
    G4LogicalVolume *ConstructDetector(G4LogicalVolume *parentVolume);
    G4LogicalVolume *ConstructEnvironment(G4LogicalVolume *parentVolume);

  private:
    bool fCheckOverlaps;
  };

} // namespace MuGeoBuilder

#endif
