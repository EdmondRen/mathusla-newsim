// G4 material
#include "G4Material.hh"
#include "G4NistManager.hh"
// G4 geometry
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
// G4 visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
// G4 constants/units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// User include:
#include "geometry/TestStand_UofT.hh"

namespace MuGeoBuilder
{
  // Use a namespace to hold all Geometry parameters
  namespace uoftdims
  {
    // 0: bar
    double bar_lenx = 100 * cm;
    double bar_leny = 4 * cm;
    double bar_lenz = 1 * cm;
    // 1: layer
    int layer_Nbars_x = 1;
    int layer_Nbars_y = 20;
    double layer_wallthick = 0.3 * cm;
    double layer_hbeam_width = 2 * cm;
    double layer_hbeam_thick = 0.5 * cm;
    double layer_lenx = bar_lenx * layer_Nbars_x + 20 * cm; //---Derived
    double layer_leny = bar_leny * layer_Nbars_y;           //---Derived
    double layer_lenz = 10 * cm;
    // 2: module
    int module_Nlayers = 4;
    double module_lenx = 140 * cm;
    double module_leny = 140 * cm;
    double module_lenz = 150 * cm;
    std::vector<int> module_layers_zdirection = {kZAxis, kZAxis, kZAxis, kZAxis};
    std::vector<int> module_layers_xdirection = {kXAxis, kYAxis, kXAxis, kYAxis};
    std::vector<double> module_layers_xoffset = {0, 1, 0, 1};
    std::vector<double> module_layers_yoffset = {0, 1, 0, 1};
    std::vector<double> module_layers_zoffset = {0, 1, 0, 1};
    double module_vbeam_width = 4 * cm;
    double module_vbeam_thick = 1 * cm;
    // Entire detector
    double detector_lenx = 120 * cm;
    double detector_leny = 90 * cm;
    double detector_lenz = 10 * cm;
    std::vector<int> detector_Nmodules_x = {1};
    std::vector<int> detector_Nmodules_y = {1};
    std::vector<double> detector_modules_xoffset = {1};
    std::vector<double> detector_modules_yoffset = {1};
    // Surroundings
    double env_ground_thickness = 120 * cm;
    double env_air_thickness = 90 * cm;
    double env_ceiling_thickness = 10 * cm;
    double env_floor_thickness = 10 * cm;
    // World
    double world_lenx = 10 * m;
    double world_leny = 10 * m;
    double world_lenz = 10 * m;
  }

  // Geometry Builder Class
  Uoft1_Builder::Uoft1_Builder(const std::string &detector_name) : Builder(detector_name)
  {
    G4cout << "Builder -------" << G4endl;
  }

  // Core function 1:
  // Construct physics volume, should return world PV
  G4VPhysicalVolume *Uoft1_Builder::Construct()
  {
    // DefineMaterials();
    G4cout << "Construct start ---------------------------------------" << G4endl;

    // World solid, logical volume, physical volume
    auto worldS = new G4Box("World",                                                                       // its name
                            uoftdims::world_lenx / 2, uoftdims::world_leny / 2, uoftdims::world_lenz / 2); // its size
    auto worldLV = new G4LogicalVolume(
        worldS,           // its solid
        Material::Vacuum, // its material
        "World");         // its name
    worldPV = new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        worldLV,         // its logical volume
        "World",         // its name
        0,               // its mother volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    // Make detector LV
    // G4LogicalVolume *detector = ConstructDetector();
    // G4LogicalVolume *environment = ConstructEnvironment();

    //
    // Calorimeter
    //
    auto calorimeterS = new G4Box("Calorimeter",                                                                 // its name
                                  uoftdims::world_lenx / 2, uoftdims::world_leny / 2, uoftdims::world_lenz / 2); // its size

    auto calorLV = new G4LogicalVolume(
        calorimeterS,   // its solid
        Material::Air,  // its material
        "Calorimeter"); // its name

    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        calorLV,         // its logical volume
        "Calorimeter",   // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    return worldPV;
  }

  // Core function 2:
  // Set the sensitive detector for this geometry

  // Helper functions:
  // No need to define material for this detector, all needed material alreade defined in MuGeoBuilder::Material
  void Uoft1_Builder::DefineMaterials() {}

  // Setup all geometry parameters
  void Uoft1_Builder::DefineGeometry()
  {
  }

  G4LogicalVolume *Uoft1_Builder::ConstructLayer()
  {
  }

  G4LogicalVolume *Uoft1_Builder::ConstructModule()
  {
  }

  G4LogicalVolume *Uoft1_Builder::ConstructDetector()
  {
  }

  G4LogicalVolume *Uoft1_Builder::ConstructEnvironment()
  {
  }

} // namespace MuGeoBuilder
