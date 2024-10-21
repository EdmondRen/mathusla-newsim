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
    std::vector<double> module_layers_xoffset = {0, 0, 0, 0};
    std::vector<double> module_layers_yoffset = {0, 0, 0, 0};
    std::vector<double> module_layers_zoffset = {0.1, 1.1, 1.6, 2.1};
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
  Uoft1_Builder::Uoft1_Builder() : Builder()
  {
    // Setup detector name and messenger
    this->DetectorName = "uoft1";
    this->fMessenger = new G4GenericMessenger(this, "/det/" + this->DetectorName, "Detector is: UofT teststand 1");
  }

  // Core function 1:
  // Construct physics volume, should return world PV
  G4VPhysicalVolume *Uoft1_Builder::Construct()
  {
    // World solid, logical volume, physical volume
    auto worldS = new G4Box("World",                                                                       // its name
                            uoftdims::world_lenx / 2, uoftdims::world_leny / 2, uoftdims::world_lenz / 2); // its size
    auto worldLV = new G4LogicalVolume(
        worldS,           // its solid
        Material::Vacuum, // its material
        "World");         // its name
    this->worldPV = new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        worldLV,         // its logical volume
        "World",         // its name
        0,               // its mother volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    // Make detector and environment materials
    G4LogicalVolume *detector = ConstructDetector(worldLV);
    // G4LogicalVolume *environment = ConstructEnvironment();

    return worldPV;
  }

  G4LogicalVolume *Uoft1_Builder::ConstructLayer(G4LogicalVolume *layerLV)
  {
    auto bar = new G4Box("bar", uoftdims::bar_lenx / 2, uoftdims::bar_leny / 2, uoftdims::bar_lenz / 2);
    auto barLV = new G4LogicalVolume(
        bar,                           // its solid
        Material::PlasticScintillator, // its material
        "bar");                        // its name

    // Fill the bars in the layer volume
    for (int i = 0; i < uoftdims::layer_Nbars_x; i++)
    {
      for (int j = 0; j < uoftdims::layer_Nbars_y; j++)
      {
        auto barPV = new G4PVPlacement(
            0, // no rotation
            G4ThreeVector(uoftdims::bar_lenx * (i - (uoftdims::layer_Nbars_x - 1) / 2.),
                          uoftdims::bar_leny * (j - (uoftdims::layer_Nbars_y - 1) / 2.),
                          0),                // offset it by corresponding bar width and length
            barLV,                           // its logical volume
            "bar",                           // its name
            layerLV,                         // its mother  volume
            false,                           // no boolean operation
            j + i * uoftdims::layer_Nbars_y, // copy number (bar number within a layer)
            fCheckOverlaps);                 // checking overlaps

        // Add this bar to the list of sensitive detectors
        allSensitiveDetectors.push_back(barPV);
      }
    }
    return 0;
  }

  G4LogicalVolume *Uoft1_Builder::ConstructModule(G4LogicalVolume *moduleLV)
  {
    // Make a layer
    auto layerS = new G4Box("layer", uoftdims::layer_lenx / 2, uoftdims::layer_leny / 2, uoftdims::layer_lenz / 2);
    auto layerLV = new G4LogicalVolume(
        layerS,        // its solid
        Material::Air, // its material
        "layer");      // its name
    ConstructLayer(layerLV);

    // Repeat the layers in module logical volume
    for (int i = 0; i < uoftdims::module_Nlayers; i++)
    {
      // Make a rotation matrix based on given directions
      std::vector<int> e1 = {0, 0, 0}, e2 = {0, 0, 0}, e3 = {0, 0, 0};
      e1[uoftdims::module_layers_xdirection[i]] = 1;
      e3[uoftdims::module_layers_zdirection[i]] = 1;
      e2[3 - uoftdims::module_layers_xdirection[i] - uoftdims::module_layers_zdirection[i]] = 1;

      auto layerPV = new G4PVPlacement(
          G4Transform3D(G4RotationMatrix(G4ThreeVector(e1[0], e1[1], e1[2]),
                                         G4ThreeVector(e2[0], e2[1], e2[2]),
                                         G4ThreeVector(e3[0], e3[1], e3[2])), // rotation
                        G4ThreeVector(uoftdims::module_layers_xoffset[i],
                                      uoftdims::module_layers_yoffset[i],
                                      -uoftdims::module_lenz/2 + uoftdims::module_layers_zoffset[i])), // offset
          layerLV,                                                          // its logical volume
          "layer",                                                          // its name
          moduleLV,                                                         // its mother volume
          false,                                                            // no boolean operation
          i,                                                                // copy number (layer number within a module)
          fCheckOverlaps);                                                  // checking overlaps
    }
    return 0;
  }

  G4LogicalVolume *Uoft1_Builder::ConstructDetector(G4LogicalVolume *detectorLV)
  {
    // Make a tower module
    auto moduleS = new G4Box("module", uoftdims::module_lenx / 2, uoftdims::module_leny / 2, uoftdims::module_lenz / 2);
    auto moduleLV = new G4LogicalVolume(
        moduleS,        // its solid
        Material::Air, // its material
        "module");      // its name
    ConstructModule(moduleLV);

    auto modulePV = new G4PVPlacement(
              G4Transform3D(G4RotationMatrix(), // rotation
                            G4ThreeVector(0,0,0)), // offset
              detectorLV,                                                          // its logical volume
              "layer",                                                          // its name
              moduleLV,                                                         // its mother volume
              false,                                                            // no boolean operation
              0,                                                                // copy number (layer number within a module)
              fCheckOverlaps);                                                  // checking overlaps    

  }

  G4LogicalVolume *Uoft1_Builder::ConstructEnvironment(G4LogicalVolume *envLV)
  {
  }

} // namespace MuGeoBuilder
