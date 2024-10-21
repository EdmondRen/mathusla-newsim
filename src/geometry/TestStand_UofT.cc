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
  namespace uoftdims
  {
    double bar_lenx = 100 * cm;
    double bar_leny = 4 * cm;
    double bar_lenz = 1 * cm;
    int layer_Nbars_x = 1;
    int layer_Nbars_y = 20;
    double layer_lenx = 120 * cm;
    double layer_leny = 90 * cm;
    double layer_lenz = 10 * cm;
    double layer_wallthick = 2 * cm;
    double layer_hbeam_width = 2 * cm;
    double layer_hbeam_thick = 0.5 * cm;
    int module_Nlayers = 4;
    double module_lenx = 140 * cm;
    double module_leny = 140 * cm;
    double module_lenz = 150 * cm;
    std::vector<int> module_layers_zdirection = {3, 3, 3, 3}; // 0,1,2 for x,y,z
    std::vector<int> module_layers_xdirection = {0, 1, 0, 1}; // 0,1,2 for x,y,z
    std::vector<double> module_layers_xoffset = {0, 1, 0, 1};
    std::vector<double> module_layers_yoffset = {0, 1, 0, 1};
    std::vector<double> module_layers_zoffset = {0, 1, 0, 1};
    double module_vbeam_width = 4 * cm;
    double module_vbeam_thick = 1 * cm;
    double detector_lenx = 120 * cm;
    double detector_leny = 90 * cm;
    double detector_lenz = 10 * cm;
    std::vector<int> detector_Nmodules_x = {1};
    std::vector<int> detector_Nmodules_y = {1};
    std::vector<double> detector_modules_xoffset = {1};
    std::vector<double> detector_modules_yoffset = {1};
    double world_lenx = 120 * cm;
    double world_leny = 90 * cm;
    double world_lenz = 10 * cm;
    double world_ground_thickness = 120 * cm;
    double world_air_thickness = 90 * cm;
    double world_ceiling_thickness = 10 * cm;
    double world_floor_thickness = 10 * cm;
  }

  // Geometry Builder Class
  Uoft1_Builder::Uoft1_Builder(const std::string &detector_name) : Builder(detector_name)
  {
  }

  // Core function 1:
  // Construct physics volume, should return world PV
  G4VPhysicalVolume *Uoft1_Builder::Construct()
  {
    DefineMaterials();
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

  G4VPhysicalVolume *ConstructLayer();
  G4VPhysicalVolume *ConstructModule();
  G4VPhysicalVolume *ConstructDetector();
  G4VPhysicalVolume *ConstructEnvironment();

} // namespace MuGeoBuilder

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MuDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4int nofLayers = 10;
  G4double absoThickness = 10. * mm;
  G4double gapThickness = 5. * mm;
  G4double calorSizeXY = 10. * cm;

  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = nofLayers * layerThickness;
  auto worldSizeXY = 1.2 * calorSizeXY;
  auto worldSizeZ = 1.2 * calorThickness;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");

  if (!defaultMaterial || !absorberMaterial || !gapMaterial)
  {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("MuDetectorConstruction::DefineVolumes()",
                "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS = new G4Box("World",                                           // its name
                          worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

  auto worldLV = new G4LogicalVolume(
      worldS,          // its solid
      defaultMaterial, // its material
      "World");        // its name

  auto worldPV = new G4PVPlacement(
      0,               // no rotation
      G4ThreeVector(), // at (0,0,0)
      worldLV,         // its logical volume
      "World",         // its name
      0,               // its mother  volume
      false,           // no boolean operation
      0,               // copy number
      fCheckOverlaps); // checking overlaps

  //
  // Calorimeter
  //
  auto calorimeterS = new G4Box("Calorimeter",                                         // its name
                                calorSizeXY / 2, calorSizeXY / 2, calorThickness / 2); // its size

  auto calorLV = new G4LogicalVolume(
      calorimeterS,    // its solid
      defaultMaterial, // its material
      "Calorimeter");  // its name

  new G4PVPlacement(
      0,               // no rotation
      G4ThreeVector(), // at (0,0,0)
      calorLV,         // its logical volume
      "Calorimeter",   // its name
      worldLV,         // its mother  volume
      false,           // no boolean operation
      0,               // copy number
      fCheckOverlaps); // checking overlaps

  //
  // Layer
  //
  auto layerS = new G4Box("Layer",                                               // its name
                          calorSizeXY / 2, calorSizeXY / 2, layerThickness / 2); // its size

  auto layerLV = new G4LogicalVolume(
      layerS,          // its solid
      defaultMaterial, // its material
      "Layer");        // its name

  new G4PVReplica(
      "Layer",         // its name
      layerLV,         // its logical volume
      calorLV,         // its mother
      kZAxis,          // axis of replication
      nofLayers,       // number of replica
      layerThickness); // witdth of replica

  //
  // Absorber
  //
  auto absorberS = new G4Box("Abso",                                               // its name
                             calorSizeXY / 2, calorSizeXY / 2, absoThickness / 2); // its size

  auto absorberLV = new G4LogicalVolume(
      absorberS,        // its solid
      absorberMaterial, // its material
      "Abso");          // its name

  fAbsorberPV = new G4PVPlacement(
      0,                                        // no rotation
      G4ThreeVector(0., 0., -gapThickness / 2), // its position
      absorberLV,                               // its logical volume
      "Abso",                                   // its name
      layerLV,                                  // its mother  volume
      false,                                    // no boolean operation
      0,                                        // copy number
      fCheckOverlaps);                          // checking overlaps

  //
  // Gap
  //
  auto gapS = new G4Box("Gap",                                               // its name
                        calorSizeXY / 2, calorSizeXY / 2, gapThickness / 2); // its size

  auto gapLV = new G4LogicalVolume(
      gapS,        // its solid
      gapMaterial, // its material
      "Gap");      // its name

  fGapPV = new G4PVPlacement(
      0,                                        // no rotation
      G4ThreeVector(0., 0., absoThickness / 2), // its position
      gapLV,                                    // its logical volume
      "Gap",                                    // its name
      layerLV,                                  // its mother  volume
      false,                                    // no boolean operation
      0,                                        // copy number
      fCheckOverlaps);                          // checking overlaps

  //
  // print parameters
  //
  G4cout
      << G4endl
      << "------------------------------------------------------------" << G4endl
      << "---> The calorimeter is " << nofLayers << " layers of: [ "
      << absoThickness / mm << "mm of " << absorberMaterial->GetName()
      << " + "
      << gapThickness / mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
      << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuDetectorConstruction::ConstructSDandField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
