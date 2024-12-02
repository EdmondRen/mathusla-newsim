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
    void ConstructSD(G4VSensitiveDetector *detector) override;

    // Core function 3:
    // (For digitizer) Get a unique detector ID for each bar based on the copy number
    // Here we define det_id as
    // det_id = (uint64)000...00AAABBBCCCXXXXX
    //    XXXXX:lower 5 digits are the bar copy number.  
    //    CCC are the "layer" copy number
    //    BBB are the "tower-module" copy number
    //    AAA are the "detector" copy number
    long long int CopynumberToDetectorID(std::vector<int> copy_numbers) override;

    // Core function 4:
    // (For digitizer) Make a map from detector ID to bar information dict
    BarPosition GetBarPosition(long long detector_id) override;    

    // Helper functions:
    // void DefineMaterials();
    // void DefineGeometry();
    G4LogicalVolume *ConstructLayer(G4LogicalVolume *parentVolume);
    G4LogicalVolume *ConstructModule(G4LogicalVolume *parentVolume);
    G4LogicalVolume *ConstructDetector(G4LogicalVolume *parentVolume);
    G4LogicalVolume *ConstructEnvironment(G4LogicalVolume *parentVolume);

  private:
    bool fCheckOverlaps;

    // For each detector ID, store a struct of BarPosition
    std::map<unsigned long long int, BarPosition> IDMaps_inLayer; // depth=0
    std::map<unsigned long long int, BarPosition> IDMaps_inDetector; // depth=2
    std::map<unsigned long long int, BarPosition> IDMaps_inWorld; // depth=3
    std::map<unsigned long long int, BarPosition> IDMaps_inTower; // depth=1
  };

} // namespace MuGeoBuilder

#endif
