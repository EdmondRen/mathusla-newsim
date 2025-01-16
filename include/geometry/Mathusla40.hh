#ifndef Mathusla40_Builder_h
#define Mathusla40_Builder_h 1

// STD libs
#include <any>

// G4 include
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// User include
#include "geometry/_GeoBuilder.hh"

namespace MuGeoBuilder
{

  namespace mu40dims
  {
    // GEO_DEPTH: number of depth. Here we stick to depth of 4: #detector, #tower-module, #layer and #bar
    extern size_t GEO_DEPTH;
    // 0: bar. x: along the bar, y: width, z: thickness
    extern double bar_lenx_real;
    extern double bar_leny_real;
    extern double bar_lenx;
    extern double bar_leny;
    extern double bar_lenz;
    // 1: layer
    extern int layer_Nbars_x;
    extern int layer_Nbars_y;
    extern double layer_wallthick;
    extern double layer_hbeam_width;
    extern double layer_hbeam_thick;
    extern double layer_lenx; //---Derived
    extern double layer_leny;           //---Derived
    extern double layer_lenz;
    // 2: tower module
    extern int module_Nlayers;
    extern double module_lenx;
    extern double module_leny;
    extern double module_lenz;
    extern std::vector<int> module_layers_zdirection;
    extern std::vector<int> module_layers_xdirection;
    extern std::vector<double> module_layers_xoffset;
    extern std::vector<double> module_layers_yoffset;
    extern std::vector<double> module_layers_zoffset;
    extern double module_vbeam_width;
    extern double module_vbeam_thick;
    // 3. Entire detector
    extern int detector_Ntowers_x;
    extern int detector_Ntowers_y;
    extern double detector_lenx;
    extern double detector_leny;
    extern double detector_lenz;
    extern double detector_decay_vol_height;
    extern std::vector<double> detector_modules_xoffset;
    extern std::vector<double> detector_modules_yoffset;
    extern std::vector<double> detector_modules_zoffset;
    extern std::vector<double> detector_ground_offset; // x and y offsets are relative to detector box center, z offset is from the bottom.
    // Surroundings
    extern double env_earth_depth_top;
    extern double env_earth_depth_mid;
    extern double env_air_depth;
    extern double env_ceiling_lenx;
    extern double env_ceiling_leny;
    extern double env_ceiling_lenz;
    extern double env_ceiling_concrete_thickness;
    extern double env_floor_iron_thickness;
    // World
    extern double world_lenx;
    extern double world_leny;
    extern double world_lenz;


    // Veto layers
    // vb: veto, bottom
    // vl: veto, left
    extern int vblayer_Nbars_x;
    extern int vblayer_Nbars_y;
    extern int vbmodule_Nlayers;
    extern std::vector<int> vbmodule_layers_zdirection;
    extern std::vector<int> vbmodule_layers_xdirection;
    extern std::vector<double> vbmodule_layers_xoffset;
    extern std::vector<double> vbmodule_layers_yoffset;
    extern std::vector<double> vbmodule_layers_zoffset;

    extern int vllayer_Nbars_x;
    extern int vllayer_Nbars_y;
    extern int vlmodule_Nlayers;
    extern std::vector<int> vlmodule_layers_zdirection;
    extern std::vector<int> vlmodule_layers_xdirection;
    extern std::vector<double> vlmodule_layers_xoffset;
    extern std::vector<double> vlmodule_layers_yoffset;
    extern std::vector<double> vlmodule_layers_zoffset;

  }

  // Geometry Builder Class
  // This class takes care of detector construction
  // as well as assigning sensitive detector to correct volumes
  class Mathusla40_Builder : public Builder
  {
  public:
    Mathusla40_Builder();

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
    long long int GetDetectorID(std::vector<int> copy_numbers, G4ThreeVector local_coord) override;

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

    void ConstructVeto(G4LogicalVolume *_worldLV);


  private:
    bool fCheckOverlaps;

    // For each detector ID, store a struct of BarPosition
    std::map<unsigned long long int, BarPosition> IDMaps_inLayer;    // depth=0
    std::map<unsigned long long int, BarPosition> IDMaps_inDetector; // depth=2
    std::map<unsigned long long int, BarPosition> IDMaps_inWorld;    // depth=3
    std::map<unsigned long long int, BarPosition> IDMaps_inTower;    // depth=1
  };

} // namespace MuGeoBuilder

#endif
