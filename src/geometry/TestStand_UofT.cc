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
#include "G4SubtractionSolid.hh"
#include "G4UserLimits.hh"
// #include "Transform3D.hh"

// G4 visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
// G4 constants/units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// User include:
#include "geometry/TestStand_UofT.hh"
#include "util.hh"

namespace MuGeoBuilder
{
  // Use a namespace to hold all Geometry parameters
  namespace uoftdims
  {
    size_t GEO_DEPTH = 4;

    // 0: bar. x: along the bar, y: width, z: thickness
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
    // 2: tower module
    int module_Nlayers = 4;
    double module_lenx = 140 * cm;
    double module_leny = 140 * cm;
    double module_lenz = 250 * cm;
    std::vector<int> module_layers_zdirection = {kZAxis, kZAxis, kZAxis, kZAxis};
    std::vector<int> module_layers_xdirection = {kXAxis, kYAxis, kXAxis, kYAxis};
    std::vector<double> module_layers_xoffset = {0, 0, 0, 0};                         // From the module x-center
    std::vector<double> module_layers_yoffset = {0, 0, 0, 0};                         // From the module y-center
    std::vector<double> module_layers_zoffset = {0.1 * m, 1.1 * m, 1.6 * m, 2.1 * m}; // From the module z-BOTTOM
    double module_vbeam_width = 4 * cm;
    double module_vbeam_thick = 1 * cm;
    double module_gap = 1 * m;
    // 3. Entire detector
    int detector_Ntowers_x = 1;
    int detector_Ntowers_y = 1;
    double detector_lenx = 140 * cm;
    double detector_leny = 140 * cm;
    double detector_lenz = 250 * cm;
    std::vector<double> detector_modules_xoffset = {0 * m};
    std::vector<double> detector_modules_yoffset = detector_modules_xoffset;
    std::vector<double> detector_ground_offset = {0, 0, 0}; // x and y offsets are relative to detector box center, z offset is from the bottom.
    // Surroundings
    double env_earth_depth_top = 0.1 * m;
    double env_earth_depth_mid = 40 * m;
    double env_air_depth = 10 * m;
    double env_ceiling_lenx = detector_lenx + 5 * m;
    double env_ceiling_leny = detector_lenx + 5 * m;
    double env_ceiling_lenz = detector_lenz + 3 * m;
    double env_ceiling_concrete_thickness = 5 * cm;
    double env_floor_iron_thickness = 2 * cm;
    // World
    double world_lenx = 100 * m;
    double world_leny = 100 * m;
    double world_lenz = 100 * m;
  }

  // Calculate detector ID based on the four copy numbers
  long long int Uoft1_Builder::GetDetectorID(std::vector<int> copy_numbers, G4ThreeVector local_coord)
  {
    (void)local_coord;
    if (copy_numbers.size() != uoftdims::GEO_DEPTH)
    {
      G4cout << " [ERROR] Geometry: The geometry depth is wrong. Please check geometry implementation." << G4endl;
      exit(0);
    }

    long long int det_id = copy_numbers[0];
    for (size_t i = 1; i < copy_numbers.size(); i++)
    {
      det_id += copy_numbers[i] * std::pow(10, 5 + i * 3); // Depth 0 takes 5 digits, the rest takes 3 digits each.
    }
    return det_id;
  }

  // Get the bar position given a detector ID
  BarPosition Uoft1_Builder::GetBarPosition(long long detector_id)
  {
    // for (auto const &[key, val] : IDMaps_inWorld)
    // print(key, val.y_side_direction);
    // print("Now asking for ID", detector_id);

    // Use find to check if the key exists
    auto it = this->IDMaps_inWorld.find(detector_id);
    if (it != this->IDMaps_inWorld.end())
    {
      return it->second;
    }
    else
    {
      print("  Error finding detector id of", detector_id);
      return BarPosition({0,0,0},{0,0,0},{0,0,0});
    }
  }

  // Return the entire map
  BarPositionMap Uoft1_Builder::GetBarPositionMap()
  {
    return this->IDMaps_inWorld;
  }

  // Geometry Builder Class
  Uoft1_Builder::Uoft1_Builder() : Builder()
  {
    // Setup detector name and messenger
    this->DetectorName = "uoft1";
    this->fMessenger = new G4GenericMessenger(this, "/det/" + this->DetectorName + "/", "Detector is: UofT teststand 1");

    // Setup the storage of bar information
    // IDMaps_inTower = new std::map<unsigned long long int, BarPosition>;
  }

  // Core function 1:
  // Construct physics volume, should return world PV
  G4VPhysicalVolume *Uoft1_Builder::Construct()
  {
    // World solid, logical volume, physical volume
    auto worldS = new G4Box("World",                                                                       // its name
                            uoftdims::world_lenx / 2, uoftdims::world_leny / 2, uoftdims::world_lenz / 2); // its size
    this->worldLV = new G4LogicalVolume(
        worldS,        // its solid
        Material::Air, // its material
        "World");      // its name
    this->worldLV->SetVisAttributes(Vis::styles["TransparentBlue"]);
    this->worldPV = new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        this->worldLV,   // its logical volume
        "World",         // its name
        0,               // its mother volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps
    (void)worldPV;

    // Make detector and environment materials, and place them in world
    ConstructDetector(this->worldLV);
    ConstructEnvironment(this->worldLV);

    return this->worldPV;
  }

  void Uoft1_Builder::ConstructSD(G4VSensitiveDetector *detector)
  {
    // Keep a copy of the pointer
    this->fdetector = detector;

    // Assign each physical volume the sensitive detector
    for (auto &bar : allSensitiveDetectors)
    {
      bar->GetLogicalVolume()->SetSensitiveDetector(detector);
    }
  }

  G4LogicalVolume *Uoft1_Builder::ConstructLayer(G4LogicalVolume *layerLV)
  {
    auto bar = new G4Box("bar", uoftdims::bar_lenx / 2, uoftdims::bar_leny / 2, uoftdims::bar_lenz / 2);
    auto barLV = new G4LogicalVolume(
        bar,                           // its solid
        Material::PlasticScintillator, // its material
        "bar");                        // its name
    barLV->SetVisAttributes(Vis::styles["SensitiveAttributes_border"]);

    // Fill the bars in the layer volume
    for (int i = 0; i < uoftdims::layer_Nbars_x; i++)
    {
      for (int j = 0; j < uoftdims::layer_Nbars_y; j++)
      {
        int bar_copy_number = j + i * uoftdims::layer_Nbars_y;

        float barcenter_x_offset = uoftdims::bar_lenx * (i - (uoftdims::layer_Nbars_x - 1) / 2.);
        float barcenter_y_offset = uoftdims::bar_leny * (j - (uoftdims::layer_Nbars_y - 1) / 2.);

        auto barPV = new G4PVPlacement(
            0, // no rotation
            G4ThreeVector(barcenter_x_offset,
                          barcenter_y_offset,
                          0), // offset it by corresponding bar width and length
            barLV,            // its logical volume
            "bar",            // its name
            layerLV,          // its mother  volume
            false,            // no boolean operation
            bar_copy_number,  // copy number (bar number within a layer)
            fCheckOverlaps);  // checking overlaps

        // Add this bar to the list of sensitive detectors
        allSensitiveDetectors.push_back(barPV);

        // Add this bar to the detector position map
        G4ThreeVector y_side_direction = G4ThreeVector(0, 1, 0);
        G4ThreeVector z_side_direction = G4ThreeVector(0, 0, 1);
        G4ThreeVector bar_center_coord = G4ThreeVector(barcenter_x_offset, barcenter_y_offset, 0);
        IDMaps_inLayer.insert({bar_copy_number, BarPosition(y_side_direction, z_side_direction, bar_center_coord)});
      }
    }

    // Add the aluminum case
    // For simplicity, just add a top and bottom plate instead of a hollow box
    auto al_case = new G4Box("al_case", uoftdims::layer_lenx / 2, uoftdims::layer_leny / 2, uoftdims::layer_wallthick / 2);
    auto al_caseLV = new G4LogicalVolume(
        al_case,            // its solid
        Material::Aluminum, // its material
        "al_case");         // its name
    new G4PVPlacement(
        0, // no rotation
        G4ThreeVector(0,
                      0,
                      uoftdims::bar_lenz / 2 + uoftdims::layer_wallthick / 2), // offset it by corresponding bar width and length
        al_caseLV,                                                             // its logical volume
        "al_case",                                                             // its name
        layerLV,                                                               // its mother  volume
        false,                                                                 // no boolean operation
        0,                                                                     // copy number
        fCheckOverlaps);                                                       // checking overlaps
    new G4PVPlacement(
        0, // no rotation
        G4ThreeVector(0,
                      0,
                      -(uoftdims::bar_lenz / 2 + uoftdims::layer_wallthick / 2)), // offset it by corresponding bar width and length
        al_caseLV,                                                                // its logical volume
        "al_case",                                                                // its name
        layerLV,                                                                  // its mother  volume
        false,                                                                    // no boolean operation
        1,                                                                        // copy number
        fCheckOverlaps);                                                          // checking overlaps

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

      auto transform_rotate = G4RotationMatrix(G4ThreeVector(e1[0], e1[1], e1[2]),
                                               G4ThreeVector(e2[0], e2[1], e2[2]),
                                               G4ThreeVector(e3[0], e3[1], e3[2]));
      auto transform_shift = G4ThreeVector(uoftdims::module_layers_xoffset[i],
                                           uoftdims::module_layers_yoffset[i],
                                           -uoftdims::module_lenz / 2 + uoftdims::module_layers_zoffset[i]); // offset

      auto layerPV = new G4PVPlacement(
          G4Transform3D(transform_rotate, transform_shift), // offset
          layerLV,                                          // its logical volume
          "layer",                                          // its name
          moduleLV,                                         // its mother volume
          false,                                            // no boolean operation
          i,                                                // copy number (layer number within a module)
          fCheckOverlaps);                                  // checking overlaps
      (void)layerPV;

      // Add the bars in this layer to the detector position map
      for (auto const &[key, val] : IDMaps_inLayer)
      {
        G4ThreeVector y_side_direction = transform_rotate * val.y_side_direction;
        G4ThreeVector z_side_direction = transform_rotate * val.z_side_direction;
        G4ThreeVector bar_center_coord = transform_rotate * val.bar_center_coord + transform_shift;
        IDMaps_inTower.insert({key + i * 1e8, BarPosition(y_side_direction, z_side_direction, bar_center_coord)});
      }
    }
    return 0;
  }

  G4LogicalVolume *Uoft1_Builder::ConstructDetector(G4LogicalVolume *_worldLV)
  {
    // Make detector logic volume
    auto detector = new G4Box("detector", uoftdims::detector_lenx / 2, uoftdims::detector_leny / 2, uoftdims::detector_lenz / 2);
    this->detectorLV = new G4LogicalVolume(
        detector,      // its solid
        Material::Air, // its material
        "detector");   // its name
    // this->detectorLV->SetVisAttributes(Vis::styles["CasingAttributes"]);

    // Place components in detector volume
    // Make a tower module
    auto moduleS = new G4Box("module", uoftdims::module_lenx / 2, uoftdims::module_leny / 2, uoftdims::module_lenz / 2);
    auto moduleLV = new G4LogicalVolume(
        moduleS,       // its solid
        Material::Air, // its material
        "module");     // its name
    ConstructModule(moduleLV);
    moduleLV->SetVisAttributes(Vis::styles["CasingAttributes"]);
    // Place all tower modules into detector
    for (int i = 0; i < uoftdims::detector_Ntowers_x; i++)
    {
      for (int j = 0; j < uoftdims::detector_Ntowers_y; j++)
      {
        int tower_copy_number = j + i * uoftdims::detector_Ntowers_y;
        auto modulePV = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(), // rotation
                          G4ThreeVector(uoftdims::detector_modules_xoffset[i],
                                        uoftdims::detector_modules_yoffset[j],
                                        0)), // offset
            moduleLV,                        // its logical volume
            "Tower module",                  // its name
            this->detectorLV,                // its mother volume
            false,                           // no boolean operation
            tower_copy_number,               // copy number (tower number within the full detector)
            fCheckOverlaps);                 // checking overlaps
        (void)modulePV;

        // Add the modules to the detector position map
        for (auto const &[key, val] : IDMaps_inTower)
        {
          G4ThreeVector bar_center_coord = val.bar_center_coord + G4ThreeVector(uoftdims::detector_modules_xoffset[i],
                                                                                uoftdims::detector_modules_yoffset[j],
                                                                                0);
          IDMaps_inDetector.insert({key + tower_copy_number * 1e8 * 1e3, BarPosition(val.y_side_direction, val.z_side_direction, bar_center_coord)});
        }
      }
    }

    // Place detector in world
    int detector_copy_number = 0;
    auto offset = G4ThreeVector(uoftdims::detector_ground_offset[0],
                                uoftdims::detector_ground_offset[1],
                                0.5 * uoftdims::detector_lenz + uoftdims::detector_ground_offset[2]);
    auto transform = G4Transform3D(G4RotationMatrix(), // rotation
                                   offset);            // offset
    auto detectorPV = new G4PVPlacement(
        transform,
        this->detectorLV,     // its logical volume
        "layer",              // its name
        _worldLV,             // its mother volume
        false,                // no boolean operation
        detector_copy_number, // copy number (detector number within the world)
        fCheckOverlaps);      // checking overlaps
    (void)detectorPV;

    // Add the detector to the detector position map
    for (auto const &[key, val] : IDMaps_inDetector)
    {
      G4ThreeVector bar_center_coord = val.bar_center_coord + offset;
      IDMaps_inWorld.insert({key + detector_copy_number * 1e8 * 1e3 * 1e3, BarPosition(val.y_side_direction, val.z_side_direction, bar_center_coord)});
    }
    return 0;
  }

  G4LogicalVolume *Uoft1_Builder::ConstructEnvironment(G4LogicalVolume *_worldLV)
  {
    // Make the detector box again, so that we can subtract it from earth/air
    auto detector = new G4Box("detector", uoftdims::detector_lenx / 2, uoftdims::detector_leny / 2, uoftdims::detector_lenz / 2);

    //
    // Earth - top
    auto earth = new G4Box("earth", uoftdims::world_lenx / 2, uoftdims::world_leny / 2, uoftdims::env_earth_depth_top / 2);
    // Subtract detector from earth (in case we want to excavate a little bit)
    auto earth_excavated =
        new G4SubtractionSolid("earth_excavated",
                               earth,
                               detector,
                               G4Transform3D(G4RotationMatrix(), // rotation
                                             G4ThreeVector(uoftdims::detector_ground_offset[0],
                                                           uoftdims::detector_ground_offset[1],
                                                           0.5 * uoftdims::env_earth_depth_top + 0.5 * uoftdims::detector_lenz + uoftdims::detector_ground_offset[2])));
    G4LogicalVolume *earthLV = new G4LogicalVolume(
        earth_excavated,     // its solid
        Material::GroundMix, // its material
        "earth");            // its name
    earthLV->SetVisAttributes(Vis::styles["TransparentBrown"]);

    // Place earth in world
    auto earthPV = new G4PVPlacement(
        G4Transform3D(G4RotationMatrix(), // rotation
                      G4ThreeVector(0, 0,
                                    -0.5 * uoftdims::env_earth_depth_top)), // offset
        earthLV,                                                            // its logical volume
        "earth",                                                            // its name
        _worldLV,                                                           // its mother volume
        false,                                                              // no boolean operation
        0,                                                                  // copy number (layer number within a module)
        fCheckOverlaps);                                                    // checking overlaps
    (void)earthPV;

    // earth - mid
    auto earth_mid = new G4Box("earth", uoftdims::world_lenx / 2, uoftdims::world_leny / 2, uoftdims::env_earth_depth_mid / 2);
    // Subtract detector from earth (in case we want to excavate a little bit)
    auto earth_excavated_mid =
        new G4SubtractionSolid("earth_excavated",
                               earth_mid,
                               detector,
                               G4Transform3D(G4RotationMatrix(), // rotation
                                             G4ThreeVector(uoftdims::detector_ground_offset[0],
                                                           uoftdims::detector_ground_offset[1],
                                                           uoftdims::env_earth_depth_top + 0.5 * uoftdims::env_earth_depth_mid + 0.5 * uoftdims::detector_lenz + uoftdims::detector_ground_offset[2])));
    G4LogicalVolume *earthLV_mid = new G4LogicalVolume(
        earth_excavated_mid, // its solid
        Material::GroundMix, // its material
        "earth");            // its name
    earthLV_mid->SetVisAttributes(Vis::styles["TransparentBrown"]);

    // Place earth in world
    auto earthPV_mid = new G4PVPlacement(
        G4Transform3D(G4RotationMatrix(), // rotation
                      G4ThreeVector(0, 0,
                                    -uoftdims::env_earth_depth_top - 0.5 * uoftdims::env_earth_depth_mid)), // offset
        earthLV_mid,                                                                                        // its logical volume
        "earth",                                                                                            // its name
        _worldLV,                                                                                           // its mother volume
        false,                                                                                              // no boolean operation
        0,                                                                                                  // copy number (layer number within a module)
        fCheckOverlaps);                                                                                    // checking overlaps
    (void)earthPV_mid;

    // Limit the step in earth mid logical volume
    G4double minEkin = 1 * MeV;   // min kinetic energy (only for charged particles)
    G4double minRange = 100 * mm; // min remaining range (only for charged particles)
    auto fStepLimit = new G4UserLimits();
    fStepLimit->SetUserMinEkine(minEkin);
    fStepLimit->SetUserMinRange(minRange);
    earthLV_mid->SetUserLimits(fStepLimit);

    //
    // ceiling
    auto ceiling_out = new G4Box("ceiling_out", uoftdims::env_ceiling_lenx / 2, uoftdims::env_ceiling_leny / 2, uoftdims::env_ceiling_lenz / 2);
    auto ceiling_inside = new G4Box("ceiling_inside", uoftdims::env_ceiling_lenx / 2 - uoftdims::env_ceiling_concrete_thickness,
                                    uoftdims::env_ceiling_leny / 2 - uoftdims::env_ceiling_concrete_thickness,
                                    uoftdims::env_ceiling_lenz / 2 - uoftdims::env_ceiling_concrete_thickness / 2);
    // Subtract detector from earth (in case we want to excavate a little bit)
    auto ceiling =
        new G4SubtractionSolid("ceiling",
                               ceiling_out,
                               ceiling_inside,
                               G4Transform3D(G4RotationMatrix(), // rotation
                                             G4ThreeVector(0, 0, -1. * uoftdims::env_ceiling_concrete_thickness)));
    G4LogicalVolume *ceilingLV = new G4LogicalVolume(
        ceiling,            // its solid
        Material::Concrete, // its material
        "ceiling");         // its name
    ceilingLV->SetVisAttributes(Vis::styles["TransparentBrown"]);

    // Place ceiling in world
    auto ceilingPV = new G4PVPlacement(
        G4Transform3D(G4RotationMatrix(), // rotation
                      G4ThreeVector(0, 0,
                                    0.5 * uoftdims::env_ceiling_lenz)), // offset
        ceilingLV,                                                      // its logical volume
        "ceiling",                                                      // its name
        _worldLV,                                                       // its mother volume
        false,                                                          // no boolean operation
        0,                                                              // copy number (layer number within a module)
        fCheckOverlaps);                                                // checking overlaps
    (void)ceilingPV;
    return 0;
  }

} // namespace MuGeoBuilder
