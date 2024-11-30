// // G4 material
// #include "G4Material.hh"
// #include "G4NistManager.hh"
// // G4 geometry
// #include "G4Box.hh"
// #include "G4LogicalVolume.hh"
// #include "G4PVPlacement.hh"
// #include "G4PVReplica.hh"
// #include "G4GlobalMagFieldMessenger.hh"
// #include "G4AutoDelete.hh"
// #include "G4GeometryManager.hh"
// #include "G4PhysicalVolumeStore.hh"
// #include "G4LogicalVolumeStore.hh"
// #include "G4SolidStore.hh"
// #include <G4SubtractionSolid.hh>

// // G4 visualization
// #include "G4VisAttributes.hh"
// #include "G4Colour.hh"
// // G4 constants/units
// #include "G4PhysicalConstants.hh"
// #include "G4SystemOfUnits.hh"

// // User include:
// #include "geometry/Mathusla40.hh"

// namespace MuGeoBuilder
// {
//   // Use a namespace to hold all Geometry parameters
//   namespace ma40dims
//   {
//     // 0: bar
//     double bar_lenx = 100 * cm;
//     double bar_leny = 4 * cm;
//     double bar_lenz = 1 * cm;
//     // 1: layer
//     int layer_Nbars_x = 1;
//     int layer_Nbars_y = 20;
//     double layer_wallthick = 0.3 * cm;
//     double layer_hbeam_width = 2 * cm;
//     double layer_hbeam_thick = 0.5 * cm;
//     double layer_lenx = bar_lenx * layer_Nbars_x + 20 * cm; //---Derived
//     double layer_leny = bar_leny * layer_Nbars_y;           //---Derived
//     double layer_lenz = 10 * cm;
//     // 2: module
//     int module_Nlayers = 4;
//     double module_lenx = 140 * cm;
//     double module_leny = 140 * cm;
//     double module_lenz = 250 * cm;
//     std::vector<int> module_layers_zdirection = {kZAxis, kZAxis, kZAxis, kZAxis};
//     std::vector<int> module_layers_xdirection = {kXAxis, kYAxis, kXAxis, kYAxis};
//     std::vector<double> module_layers_xoffset = {0, 0, 0, 0};                         // From the module x-center
//     std::vector<double> module_layers_yoffset = {0, 0, 0, 0};                         // From the module y-center
//     std::vector<double> module_layers_zoffset = {0.1 * m, 1.1 * m, 1.6 * m, 2.1 * m}; // From the module z-BOTTOM
//     double module_vbeam_width = 4 * cm;
//     double module_vbeam_thick = 1 * cm;
//     double module_gap = 1 * m;
//     // Entire detector
//     int detector_Ntowers_x = 1;
//     int detector_Ntowers_y = 1;
//     double detector_lenx = 140 * cm;
//     double detector_leny = 140 * cm;
//     double detector_lenz = 250 * cm;
//     std::vector<double> detector_modules_xoffset = {0 * m};
//     std::vector<double> detector_modules_yoffset = detector_modules_xoffset;
//     std::vector<double> detector_ground_offset = {0, 0, 0}; // x and y offsets are relative to detector box center, z offset is from the bottom.
//     // Surroundings
//     double env_earth_depth = 10 * m;
//     double env_air_depth = 10 * m;
//     double env_ceiling_lenx = detector_lenx + 5 * m;
//     double env_ceiling_leny = detector_lenx + 5 * m;
//     double env_ceiling_lenz = detector_lenz + 3 * m;
//     double env_ceiling_concrete_thickness = 40 * cm;
//     double env_floor_iron_thickness = 2 * cm;
//     // World
//     double world_lenx = 10 * m;
//     double world_leny = 10 * m;
//     double world_lenz = 40 * m;
//   }

//   // Geometry Builder Class
//   Mathusla40_Builder::Mathusla40_Builder() : Builder()
//   {
//     // Setup detector name and messenger
//     this->DetectorName = "mathusla40";
//     this->fMessenger = new G4GenericMessenger(this, "/det/" + this->DetectorName + "/", "Detector is: UofT teststand 1");
//   }

//   // Core function 1:
//   // Construct physics volume, should return world PV
//   G4VPhysicalVolume *Mathusla40_Builder::Construct()
//   {
//     // World solid, logical volume, physical volume
//     auto worldS = new G4Box("World",                                                                       // its name
//                             ma40dims::world_lenx / 2, ma40dims::world_leny / 2, ma40dims::world_lenz / 2); // its size
//     this->worldLV = new G4LogicalVolume(
//         worldS,           // its solid
//         Material::Air, // its material
//         "World");         // its name
//     this->worldLV->SetVisAttributes(Vis::styles["TransparentBlue"]);
//     this->worldPV = new G4PVPlacement(
//         0,               // no rotation
//         G4ThreeVector(), // at (0,0,0)
//         this->worldLV,   // its logical volume
//         "World",         // its name
//         0,               // its mother volume
//         false,           // no boolean operation
//         0,               // copy number
//         fCheckOverlaps); // checking overlaps

//     // Make detector and environment materials, and place them in world
//     ConstructDetector(this->worldLV);
//     ConstructEnvironment(this->worldLV);

//     return this->worldPV;
//   }

//   void Mathusla40_Builder::ConstructSD(G4VSensitiveDetector *detector) 
//   { 
//     // Keep a copy of the pointer
//     this->fdetector = detector;

//     // Assign each physical volume the sensitive detector
//     for (auto& bar : allSensitiveDetectors)
//     {
//       bar->GetLogicalVolume()->SetSensitiveDetector(detector);
//     }
//   }

//   G4LogicalVolume *Mathusla40_Builder::ConstructLayer(G4LogicalVolume *layerLV)
//   {
//     auto bar = new G4Box("bar", ma40dims::bar_lenx / 2, ma40dims::bar_leny / 2, ma40dims::bar_lenz / 2);
//     auto barLV = new G4LogicalVolume(
//         bar,                           // its solid
//         Material::PlasticScintillator, // its material
//         "bar");                        // its name
//     barLV->SetVisAttributes(Vis::styles["SensitiveAttributes_border"]);

//     // Fill the bars in the layer volume
//     for (int i = 0; i < ma40dims::layer_Nbars_x; i++)
//     {
//       for (int j = 0; j < ma40dims::layer_Nbars_y; j++)
//       {
//         auto barPV = new G4PVPlacement(
//             0, // no rotation
//             G4ThreeVector(ma40dims::bar_lenx * (i - (ma40dims::layer_Nbars_x - 1) / 2.),
//                           ma40dims::bar_leny * (j - (ma40dims::layer_Nbars_y - 1) / 2.),
//                           0),                // offset it by corresponding bar width and length
//             barLV,                           // its logical volume
//             "bar",                           // its name
//             layerLV,                         // its mother  volume
//             false,                           // no boolean operation
//             j + i * ma40dims::layer_Nbars_y, // copy number (bar number within a layer)
//             fCheckOverlaps);                 // checking overlaps

//         // Add this bar to the list of sensitive detectors
//         allSensitiveDetectors.push_back(barPV);
//       }
//     }
//     return 0;
//   }

//   G4LogicalVolume *Mathusla40_Builder::ConstructModule(G4LogicalVolume *moduleLV)
//   {
//     // Make a layer
//     auto layerS = new G4Box("layer", ma40dims::layer_lenx / 2, ma40dims::layer_leny / 2, ma40dims::layer_lenz / 2);
//     auto layerLV = new G4LogicalVolume(
//         layerS,        // its solid
//         Material::Air, // its material
//         "layer");      // its name
//     ConstructLayer(layerLV);

//     // Repeat the layers in module logical volume
//     for (int i = 0; i < ma40dims::module_Nlayers; i++)
//     {
//       // Make a rotation matrix based on given directions
//       std::vector<int> e1 = {0, 0, 0}, e2 = {0, 0, 0}, e3 = {0, 0, 0};
//       e1[ma40dims::module_layers_xdirection[i]] = 1;
//       e3[ma40dims::module_layers_zdirection[i]] = 1;
//       e2[3 - ma40dims::module_layers_xdirection[i] - ma40dims::module_layers_zdirection[i]] = 1;

//       auto layerPV = new G4PVPlacement(
//           G4Transform3D(G4RotationMatrix(G4ThreeVector(e1[0], e1[1], e1[2]),
//                                          G4ThreeVector(e2[0], e2[1], e2[2]),
//                                          G4ThreeVector(e3[0], e3[1], e3[2])), // rotation
//                         G4ThreeVector(ma40dims::module_layers_xoffset[i],
//                                       ma40dims::module_layers_yoffset[i],
//                                       -ma40dims::module_lenz / 2 + ma40dims::module_layers_zoffset[i])), // offset
//           layerLV,                                                                                       // its logical volume
//           "layer",                                                                                       // its name
//           moduleLV,                                                                                      // its mother volume
//           false,                                                                                         // no boolean operation
//           i,                                                                                             // copy number (layer number within a module)
//           fCheckOverlaps);                                                                               // checking overlaps
//       (void)layerPV;
//     }
//     return 0;
//   }

//   G4LogicalVolume *Mathusla40_Builder::ConstructDetector(G4LogicalVolume *worldLV)
//   {
//     // Make detector logic volume
//     auto detector = new G4Box("detector", ma40dims::detector_lenx / 2, ma40dims::detector_leny / 2, ma40dims::detector_lenz / 2);
//     this->detectorLV = new G4LogicalVolume(
//         detector,      // its solid
//         Material::Air, // its material
//         "detector");   // its name
//     // this->detectorLV->SetVisAttributes(Vis::styles["CasingAttributes"]);

//     // Place detector in world
//     auto detectorPV = new G4PVPlacement(
//         G4Transform3D(G4RotationMatrix(), // rotation
//                       G4ThreeVector(ma40dims::detector_ground_offset[0],
//                                     ma40dims::detector_ground_offset[1],
//                                     0.5 * ma40dims::detector_lenz + ma40dims::detector_ground_offset[2])), // offset
//         this->detectorLV,                                                                                  // its logical volume
//         "layer",                                                                                           // its name
//         worldLV,                                                                                           // its mother volume
//         false,                                                                                             // no boolean operation
//         0,                                                                                                 // copy number (layer number within a module)
//         fCheckOverlaps);                                                                                   // checking overlaps
//     (void)detectorPV;

//     // Place components in detector volume
//     // Make a tower module
//     auto moduleS = new G4Box("module", ma40dims::module_lenx / 2, ma40dims::module_leny / 2, ma40dims::module_lenz / 2);
//     auto moduleLV = new G4LogicalVolume(
//         moduleS,       // its solid
//         Material::Air, // its material
//         "module");     // its name
//     ConstructModule(moduleLV);
//     moduleLV->SetVisAttributes(Vis::styles["CasingAttributes"]);
//     // Place all tower modules into detector
//     for (int i = 0; i < ma40dims::detector_Ntowers_x; i++)
//     {
//       for (int j = 0; j < ma40dims::detector_Ntowers_y; j++)
//       {
//         auto modulePV = new G4PVPlacement(
//             G4Transform3D(G4RotationMatrix(), // rotation
//                           G4ThreeVector(ma40dims::detector_modules_xoffset[i],
//                                         ma40dims::detector_modules_yoffset[j],
//                                         0)), // offset
//             moduleLV,                        // its logical volume
//             "layer",                         // its name
//             this->detectorLV,                // its mother volume
//             false,                           // no boolean operation
//             0,                               // copy number (layer number within a module)
//             fCheckOverlaps);                 // checking overlaps
//         (void)modulePV;
//       }
//     }

//     return 0;
//   }

//   G4LogicalVolume *Mathusla40_Builder::ConstructEnvironment(G4LogicalVolume *_worldLV)
//   {
//     // Make the detector box again, so that we can subtract it from earth/air
//     auto detector = new G4Box("detector", ma40dims::detector_lenx / 2, ma40dims::detector_leny / 2, ma40dims::detector_lenz / 2);

//     //
//     // Earth
//     auto earth = new G4Box("earth", ma40dims::world_lenx / 2, ma40dims::world_leny / 2, ma40dims::env_earth_depth / 2);
//     // Subtract detector from earth (in case we want to excavate a little bit)
//     auto earth_excavated =
//         new G4SubtractionSolid("earth_excavated",
//                                earth,
//                                detector,
//                                G4Transform3D(G4RotationMatrix(), // rotation
//                                              G4ThreeVector(ma40dims::detector_ground_offset[0],
//                                                            ma40dims::detector_ground_offset[1],
//                                                            0.5 * ma40dims::env_earth_depth + 0.5 * ma40dims::detector_lenz + ma40dims::detector_ground_offset[2])));
//     G4LogicalVolume *earthLV = new G4LogicalVolume(
//         earth_excavated,     // its solid
//         Material::GroundMix, // its material
//         "earth");            // its name
//     earthLV->SetVisAttributes(Vis::styles["TransparentBrown"]);

//     // Place earth in world
//     auto earthPV = new G4PVPlacement(
//         G4Transform3D(G4RotationMatrix(), // rotation
//                       G4ThreeVector(0, 0,
//                                     -0.5 * ma40dims::env_earth_depth)), // offset
//         earthLV,                                                        // its logical volume
//         "earth",                                                        // its name
//         _worldLV,                                                        // its mother volume
//         false,                                                          // no boolean operation
//         0,                                                              // copy number (layer number within a module)
//         fCheckOverlaps);                                                // checking overlaps

//     //
//     // ceiling
//     auto ceiling_out = new G4Box("ceiling_out", ma40dims::env_ceiling_lenx / 2, ma40dims::env_ceiling_leny / 2, ma40dims::env_ceiling_lenz / 2);
//     auto ceiling_inside = new G4Box("ceiling_inside", ma40dims::env_ceiling_lenx / 2 - ma40dims::env_ceiling_concrete_thickness,
//                                     ma40dims::env_ceiling_leny / 2 - ma40dims::env_ceiling_concrete_thickness,
//                                     ma40dims::env_ceiling_lenz / 2 - ma40dims::env_ceiling_concrete_thickness /2);
//     // Subtract detector from earth (in case we want to excavate a little bit)
//     auto ceiling =
//         new G4SubtractionSolid("ceiling",
//                                ceiling_out,
//                                ceiling_inside,
//                                G4Transform3D(G4RotationMatrix(), // rotation
//                                              G4ThreeVector(0,0, -1.*ma40dims::env_ceiling_concrete_thickness)));
//     G4LogicalVolume *ceilingLV = new G4LogicalVolume(
//         ceiling,     // its solid
//         Material::Concrete, // its material
//         "ceiling");            // its name
//     ceilingLV->SetVisAttributes(Vis::styles["TransparentBrown"]);

//     // Place ceiling in world
//     auto ceilingPV = new G4PVPlacement(
//         G4Transform3D(G4RotationMatrix(), // rotation
//                       G4ThreeVector(0, 0,
//                                     0.5 * ma40dims::env_ceiling_lenz)), // offset
//         ceilingLV,                                                        // its logical volume
//         "ceiling",                                                        // its name
//         _worldLV,                                                        // its mother volume
//         false,                                                          // no boolean operation
//         0,                                                              // copy number (layer number within a module)
//         fCheckOverlaps);                                                // checking overlaps
//   }

// } // namespace MuGeoBuilder
