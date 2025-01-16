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
#include "geometry/Mathusla40.hh"
#include "util.hh"

namespace MuGeoBuilder
{
    // Use a namespace to hold all Geometry parameters
    namespace mu40dims
    {
        size_t GEO_DEPTH = 4;

        // 0: bar. x: along the bar, y: width, z: thickness
        double bar_lenx_real = 2.25 * m;
        double bar_leny_real = 3.5 * cm;
        double bar_lenx = bar_lenx_real * floor(900 * cm / bar_lenx_real);
        double bar_leny = bar_leny_real * floor(900 * cm / bar_leny_real);
        double bar_lenz = 1 * cm;
        // 1: layer
        int layer_Nbars_x = 1;
        int layer_Nbars_y = 1;
        int layer_Nbars_x_real = round(bar_lenx / bar_lenx_real);
        int layer_Nbars_y_real = round(bar_leny / bar_leny_real);
        double layer_wallthick = 0.3 * cm;
        double layer_hbeam_width = 2 * cm;
        double layer_hbeam_thick = 0.5 * cm;
        double layer_lenx = bar_lenx * layer_Nbars_x + 20 * cm; //---Derived
        double layer_leny = bar_leny * layer_Nbars_y;           //---Derived
        double layer_lenz = 10 * cm;
        // 2: tower module
        int module_Nlayers = 6;
        double module_lenx = 10 * m;
        double module_leny = 10 * m;
        double module_lenz = 4.1 * m;
        double module_lgap = 0.8 * m; // Gap between layers
        std::vector<int> module_layers_zdirection = {kZAxis, kZAxis, kZAxis, kZAxis, kZAxis, kZAxis};
        std::vector<int> module_layers_xdirection = {kXAxis, kYAxis, kXAxis, kYAxis, kXAxis, kYAxis};
        std::vector<double> module_layers_xoffset = {0, 0, 0, 0, 0};                                                                                                                                                                            // From the module x-center
        std::vector<double> module_layers_yoffset = {0, 0, 0, 0, 0};                                                                                                                                                                            // From the module y-center
        std::vector<double> module_layers_zoffset = {layer_lenz / 2, layer_lenz / 2 + module_lgap * 1, layer_lenz / 2 + module_lgap * 2, layer_lenz / 2 + module_lgap * 3, layer_lenz / 2 + module_lgap * 4, layer_lenz / 2 + module_lgap * 5}; // From the module z-BOTTOM
        double module_vbeam_width = 10 * cm;
        double module_vbeam_thick = 1 * cm;
        // 3. Entire detector
        int detector_Ntowers_x = 4;
        int detector_Ntowers_y = 4;
        double detector_decay_vol_height = 11 * m;
        double detector_lenx = detector_Ntowers_x * module_lenx + 10 * m;
        double detector_leny = detector_Ntowers_y * module_leny;
        double detector_lenz = module_lenz;
        std::vector<double> detector_modules_xoffset = {-1.5 * module_lenx, -0.5 * module_lenx, 0.5 * module_lenx, 1.5 * module_lenx};
        std::vector<double> detector_modules_yoffset = detector_modules_xoffset;
        std::vector<double> detector_ground_offset = {0, 0, detector_decay_vol_height + module_lgap + layer_lenz}; // x and y offsets are relative to detector box center, z offset is from the bottom.
        // Surroundings
        double env_earth_depth_top = 0.1 * m;
        double env_earth_depth_mid = 100 * m;
        double env_air_depth = 20 * m;
        double env_ceiling_lenx = detector_lenx + 10 * m;
        double env_ceiling_leny = detector_lenx + 10 * m;
        double env_ceiling_lenz = detector_lenz + 3 * m + detector_ground_offset[2];
        double env_ceiling_concrete_thickness = 5 * cm;
        double env_floor_iron_thickness = 2 * cm;
        // World
        double world_lenx = 120 * m;
        double world_leny = 100 * m;
        double world_lenz = 100 * m;

        // Veto layers
        // Those are treated differently.
        // vf_: veto_floor
        // vw_: veto_wall
        double vf_panel_lenx = bar_lenx_real * floor(detector_Ntowers_x * module_lenx / bar_lenx_real);
        double vf_panel_leny = bar_leny_real * floor(detector_Ntowers_y * module_leny / bar_leny_real);
        int vf_layer_Nbars_x_real = round(vf_panel_lenx / bar_lenx_real);
        int vf_layer_Nbars_y_real = round(vf_panel_leny / bar_leny_real);
        double vf_layer_lenx = bar_lenx * layer_Nbars_x + 20 * cm;
        double vf_layer_leny = bar_leny * layer_Nbars_y;

        double veto_wall_height = detector_decay_vol_height + module_lgap + layer_lenz;
        int vw_nbars_y1 = floor(detector_leny / bar_lenx_real);
        int vw_nbars_z1 = floor(veto_wall_height / bar_leny_real);
        int vw_nbars_y2 = floor(detector_leny / bar_leny_real);
        int vw_nbars_z2 = floor(veto_wall_height / bar_lenx_real);
        double vw_panel_leny1 = bar_lenx_real * vw_nbars_y1;
        double vw_panel_lenz1 = bar_leny_real * vw_nbars_z1;
        double vw_panel_leny2 = bar_leny_real * vw_nbars_y2;
        double vw_panel_lenz2 = bar_lenx_real * vw_nbars_z2;

    }

    // Calculate detector ID based on the four copy numbers
    long long int Mathusla40_Builder::GetDetectorID(std::vector<int> copy_numbers, G4ThreeVector local_coord)
    {
        if (copy_numbers.size() != mu40dims::GEO_DEPTH)
        {
            G4cout << " [ERROR] Geometry: The geometry depth is wrong. Please check geometry implementation." << G4endl;
            exit(0);
        }

        // long long int det_id = copy_numbers[0];
        long long int det_id = 0;
        // Now we need to manually calculate which bar it is in based on the local coordinate.
        int nx = floor((local_coord.X + mu40dims::bar_lenx * 0.5) / mu40dims::bar_lenx_real);
        int ny = floor((local_coord.Y + mu40dims::bar_leny * 0.5) / mu40dims::bar_leny_real);
        int bar_copy_number = ny + nx * mu40dims::layer_Nbars_y_real;

        for (size_t i = 1; i < copy_numbers.size(); i++)
        {
            det_id += copy_numbers[i] * std::pow(10, 5 + i * 3); // Depth 0 takes 5 digits, the rest takes 3 digits each.
        }
        return det_id;
    }

    // Get the bar position given a detector ID
    BarPosition Mathusla40_Builder::GetBarPosition(long long detector_id)
    {
        // for (auto const &[key, val] : IDMaps_inWorld)
        // {
        // print(key, val.y_side_direction);
        // }

        // print("Now asking for ID", detector_id);
        return this->IDMaps_inWorld.at(detector_id);
    }

    // Geometry Builder Class
    Mathusla40_Builder::Mathusla40_Builder() : Builder()
    {
        // Setup detector name and messenger
        this->DetectorName = "mu40v0";
        this->fMessenger = new G4GenericMessenger(this, "/det/" + this->DetectorName + "/", "Detector is: Mathusla40 version 0");

        // Setup the storage of bar information
        // IDMaps_inTower = new std::map<unsigned long long int, BarPosition>;
        // IDMaps_inDetector = new std::map<unsigned long long int, BarPosition>;
        // IDMaps_inWorld = new std::map<unsigned long long int, BarPosition>;
        // IDMaps_inTower = new std::map<unsigned long long int, BarPosition>;
    }

    // Core function 1:
    // Construct physics volume, should return world PV
    G4VPhysicalVolume *Mathusla40_Builder::Construct()
    {
        // World solid, logical volume, physical volume
        auto worldS = new G4Box("World",                                                                       // its name
                                mu40dims::world_lenx / 2, mu40dims::world_leny / 2, mu40dims::world_lenz / 2); // its size
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

    void Mathusla40_Builder::ConstructSD(G4VSensitiveDetector *detector)
    {
        // Keep a copy of the pointer
        this->fdetector = detector;

        // Assign each physical volume the sensitive detector
        for (auto &bar : allSensitiveDetectors)
        {
            bar->GetLogicalVolume()->SetSensitiveDetector(detector);
        }
    }

    G4LogicalVolume *Mathusla40_Builder::ConstructLayer(G4LogicalVolume *layerLV)
    {
        auto bar = new G4Box("bar", mu40dims::bar_lenx / 2, mu40dims::bar_leny / 2, mu40dims::bar_lenz / 2);
        auto barLV = new G4LogicalVolume(
            bar,                           // its solid
            Material::PlasticScintillator, // its material
            "bar");                        // its name
        barLV->SetVisAttributes(Vis::styles["SensitiveAttributes1"]);

        // Fill the bars in the layer volume
        for (int i = 0; i < mu40dims::layer_Nbars_x; i++)
        {
            for (int j = 0; j < mu40dims::layer_Nbars_y; j++)
            {
                int bar_copy_number = j + i * mu40dims::layer_Nbars_y;

                float barcenter_x_offset = mu40dims::bar_lenx * (i - (mu40dims::layer_Nbars_x - 1) / 2.);
                float barcenter_y_offset = mu40dims::bar_leny * (j - (mu40dims::layer_Nbars_y - 1) / 2.);

                auto barPV = new G4PVPlacement(
                    0, // no rotation
                    G4ThreeVector(barcenter_x_offset,
                                  barcenter_y_offset,
                                  0), // offset it by corresponding bar width and length
                    barLV,            // its logical volume
                    "bar",            // its name
                    layerLV,          // its mother volume
                    false,            // no boolean operation
                    bar_copy_number,  // copy number (bar number within a layer)
                    fCheckOverlaps);  // checking overlaps

                // Add this bar to the list of sensitive detectors
                allSensitiveDetectors.push_back(barPV);
            }
        }

        // Need to manually calculate the real copy number of each bar
        for (int i = 0; i < mu40dims::layer_Nbars_x_real; i++)
        {
            for (int j = 0; j < mu40dims::layer_Nbars_y_real; j++)
            {
                int bar_copy_number = j + i * mu40dims::layer_Nbars_y_real;

                float barcenter_x_offset = mu40dims::bar_lenx_real * (i - (mu40dims::layer_Nbars_x_real - 1) / 2.);
                float barcenter_y_offset = mu40dims::bar_leny_real * (j - (mu40dims::layer_Nbars_y_real - 1) / 2.);
                // Add this bar to the detector position map
                G4ThreeVector y_side_direction = G4ThreeVector(0, 1, 0);
                G4ThreeVector z_side_direction = G4ThreeVector(0, 0, 1);
                G4ThreeVector bar_center_coord = G4ThreeVector(barcenter_x_offset, barcenter_y_offset, 0);
                IDMaps_inLayer.insert({bar_copy_number, BarPosition(y_side_direction, z_side_direction, bar_center_coord)});
            }
        }

        // Add the aluminum case
        // For simplicity, just add a top and bottom plate instead of a hollow box
        auto al_case = new G4Box("al_case", mu40dims::layer_lenx / 2, mu40dims::layer_leny / 2, mu40dims::layer_wallthick / 2);
        auto al_caseLV = new G4LogicalVolume(
            al_case,            // its solid
            Material::Aluminum, // its material
            "al_case");         // its name
        auto al_case1PV = new G4PVPlacement(
            0, // no rotation
            G4ThreeVector(0,
                          0,
                          mu40dims::bar_lenz / 2 + mu40dims::layer_wallthick / 2), // offset it by corresponding bar width and length
            al_caseLV,                                                             // its logical volume
            "al_case",                                                             // its name
            layerLV,                                                               // its mother  volume
            false,                                                                 // no boolean operation
            0,                                                                     // copy number
            fCheckOverlaps);                                                       // checking overlaps
        auto al_case2PV = new G4PVPlacement(
            0, // no rotation
            G4ThreeVector(0,
                          0,
                          -(mu40dims::bar_lenz / 2 + mu40dims::layer_wallthick / 2)), // offset it by corresponding bar width and length
            al_caseLV,                                                                // its logical volume
            "al_case",                                                                // its name
            layerLV,                                                                  // its mother  volume
            false,                                                                    // no boolean operation
            1,                                                                        // copy number
            fCheckOverlaps);                                                          // checking overlaps

        return 0;
    }

    G4LogicalVolume *Mathusla40_Builder::ConstructModule(G4LogicalVolume *moduleLV)
    {
        // Make a layer
        auto layerS = new G4Box("layer", mu40dims::layer_lenx / 2, mu40dims::layer_leny / 2, mu40dims::layer_lenz / 2);
        auto layerLV = new G4LogicalVolume(
            layerS,        // its solid
            Material::Air, // its material
            "layer");      // its name
        ConstructLayer(layerLV);

        // Repeat the layers in module logical volume
        for (int i = 0; i < mu40dims::module_Nlayers; i++)
        {
            // Make a rotation matrix based on given directions
            std::vector<int> e1 = {0, 0, 0}, e2 = {0, 0, 0}, e3 = {0, 0, 0};
            e1[mu40dims::module_layers_xdirection[i]] = 1;
            e3[mu40dims::module_layers_zdirection[i]] = 1;
            e2[3 - mu40dims::module_layers_xdirection[i] - mu40dims::module_layers_zdirection[i]] = 1;

            auto transform_rotate = G4RotationMatrix(G4ThreeVector(e1[0], e1[1], e1[2]),
                                                     G4ThreeVector(e2[0], e2[1], e2[2]),
                                                     G4ThreeVector(e3[0], e3[1], e3[2]));
            auto transform_shift = G4ThreeVector(mu40dims::module_layers_xoffset[i],
                                                 mu40dims::module_layers_yoffset[i],
                                                 -mu40dims::module_lenz / 2 + mu40dims::module_layers_zoffset[i]); // offset

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

    G4LogicalVolume *Mathusla40_Builder::ConstructDetector(G4LogicalVolume *_worldLV)
    {
        // Make detector logic volume
        auto detector = new G4Box("detector", mu40dims::detector_lenx / 2, mu40dims::detector_leny / 2, mu40dims::detector_lenz / 2);
        this->detectorLV = new G4LogicalVolume(
            detector,      // its solid
            Material::Air, // its material
            "detector");   // its name
        this->detectorLV->SetVisAttributes(Vis::styles["Invisible"]);

        // Place components in detector volume
        // Make a tower module
        auto moduleS = new G4Box("module", mu40dims::module_lenx / 2, mu40dims::module_leny / 2, mu40dims::module_lenz / 2);
        auto moduleLV = new G4LogicalVolume(
            moduleS,       // its solid
            Material::Air, // its material
            "module");     // its name
        ConstructModule(moduleLV);
        moduleLV->SetVisAttributes(Vis::styles["CasingAttributes"]);
        // Place all tower modules into detector
        for (int i = 0; i < mu40dims::detector_Ntowers_x; i++)
        {
            for (int j = 0; j < mu40dims::detector_Ntowers_y; j++)
            {
                int tower_copy_number = j + i * mu40dims::detector_Ntowers_y;
                auto modulePV = new G4PVPlacement(
                    G4Transform3D(G4RotationMatrix(), // rotation
                                  G4ThreeVector(mu40dims::detector_modules_xoffset[i],
                                                mu40dims::detector_modules_yoffset[j],
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
                    G4ThreeVector bar_center_coord = val.bar_center_coord + G4ThreeVector(mu40dims::detector_modules_xoffset[i],
                                                                                          mu40dims::detector_modules_yoffset[j],
                                                                                          0);
                    IDMaps_inDetector.insert({key + tower_copy_number * 1e8 * 1e3, BarPosition(val.y_side_direction, val.z_side_direction, bar_center_coord)});
                }
            }
        }
        // Backwall detectors
        for (int i = 0; i < mu40dims::detector_Ntowers_y; i++)
        {
            // Define the rotation around the Y-axis (using the angle in radians)
            G4RotationMatrix rotation;
            rotation.rotateY(90.0 * deg);

            int tower_copy_number = i + mu40dims::detector_Ntowers_x * mu40dims::detector_Ntowers_y;
            ;
            auto modulePV = new G4PVPlacement(
                G4Transform3D(rotation, // rotation
                              G4ThreeVector(mu40dims::detector_modules_xoffset.back() + mu40dims::module_lenx * 0.5 + mu40dims::module_lenz * 0.5 + 10 * cm,
                                            mu40dims::detector_modules_yoffset[i],
                                            -mu40dims::module_leny * 0.5 - mu40dims::module_lenz * 0.5)), // offset
                moduleLV,                                                                                 // its logical volume
                "Tower module",                                                                           // its name
                this->detectorLV,                                                                         // its mother volume
                false,                                                                                    // no boolean operation
                tower_copy_number,                                                                        // copy number (tower number within the full detector)
                fCheckOverlaps);                                                                          // checking overlaps
            (void)modulePV;

            // Add the modules to the detector position map
            for (auto const &[key, val] : IDMaps_inTower)
            {
                G4ThreeVector bar_center_coord = val.bar_center_coord + G4ThreeVector(mu40dims::detector_modules_xoffset.back() + mu40dims::module_lenx * 0.5 + mu40dims::module_lenz * 0.5 + 10 * cm,
                                                                                      mu40dims::detector_modules_yoffset[i],
                                                                                      0);
                IDMaps_inDetector.insert({key + tower_copy_number * 1e8 * 1e3,
                                          BarPosition(rotation * val.y_side_direction, rotation * val.z_side_direction, bar_center_coord)});
            }
        }

        // Veto layers

        ConstructVeto(_worldLV);

        // Place detector in world
        int detector_copy_number = 0;
        auto offset = G4ThreeVector(mu40dims::detector_ground_offset[0],
                                    mu40dims::detector_ground_offset[1],
                                    0.5 * mu40dims::detector_lenz + mu40dims::detector_ground_offset[2]);
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

    void Mathusla40_Builder::ConstructVeto(G4LogicalVolume *_worldLV)
    {
        // Make detector logic volume
        auto detector_z_thickness = (mu40dims::module_lgap + mu40dims::layer_lenz * 2);
        auto detector_veto_floor = new G4Box("detector", mu40dims::detector_lenx / 2, mu40dims::detector_leny / 2, detector_z_thickness / 2);
        auto detector_veto_floorLV = new G4LogicalVolume(
            detector_veto_floor, // its solid
            Material::Air,       // its material
            "detector");         // its name
        detector_veto_floorLV->SetVisAttributes(Vis::styles["Invisible"]);
        // Place detector in world
        auto detector_veto_floorPV = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(),                               // rotation
                          G4ThreeVector(0, 0, detector_z_thickness * 0.5)), // offset),
            detector_veto_floorLV,                                          // its logical volume
            "Floor veto detector",                                          // its name
            _worldLV,                                                       // its mother volume
            false,                                                          // no boolean operation
            1,                                                              // copy number (detector number within the world)
            fCheckOverlaps);                                                // checking overlaps

        // Make a module
        auto module_veto_floor = new G4Box("module", mu40dims::vf_panel_lenx / 2, mu40dims::vf_panel_leny / 2, detector_z_thickness / 2);
        auto module_veto_floorLV = new G4LogicalVolume(
            module_veto_floor, // its solid
            Material::Air,     // its material
            "module");         // its name
        module_veto_floorLV->SetVisAttributes(Vis::styles["CasingAttributes"]);
        // Place module in detector
        auto module_veto_floorPV = new G4PVPlacement(
            G4Transform3D(),
            module_veto_floorLV,   // its logical volume
            "Floor veto module",   // its name
            detector_veto_floorLV, // its mother volume
            false,                 // no boolean operation
            0,                     // copy number (detector number within the world)
            fCheckOverlaps);       // checking overlaps

        // Put layers into module
        // Make a layer
        auto layer_veto_floor = new G4Box("layer", mu40dims::vf_panel_lenx / 2, mu40dims::vf_panel_leny / 2, mu40dims::layer_lenz / 2);
        auto layer_veto_floorLV = new G4LogicalVolume(
            layer_veto_floor, // its solid
            Material::Air,    // its material
            "layer");         // its name

        // Make a scintillator
        auto bar_veto_floor = new G4Box("bar", mu40dims::vf_panel_lenx / 2, mu40dims::vf_panel_leny / 2, mu40dims::bar_lenz / 2);
        auto bar_veto_floorLV = new G4LogicalVolume(
            bar_veto_floor,                // its solid
            Material::PlasticScintillator, // its material
            "bar");                        // its name
        bar_veto_floorLV->SetVisAttributes(Vis::styles["SensitiveAttributes2"]);
        // Put scintillators in layer
        auto bar_veto_floor1PV = new G4PVPlacement(
            G4Transform3D(),
            bar_veto_floorLV,   // its logical volume
            "Bar",              // its name
            layer_veto_floorLV, // its mother volume
            false,              // no boolean operation
            0,                  // copy number (detector number within the world)
            fCheckOverlaps);    // checking overlaps


        // Add the aluminum case
        // For simplicity, just add a top and bottom plate instead of a hollow box
        auto al_case = new G4Box("al_case", mu40dims::vf_panel_lenx / 2, mu40dims::vf_panel_leny / 2, mu40dims::layer_wallthick / 2);
        auto al_caseLV = new G4LogicalVolume(
            al_case,            // its solid
            Material::Aluminum, // its material
            "al_case");         // its name
        auto al_case1PV = new G4PVPlacement(
            0, // no rotation
            G4ThreeVector(0,
                          0,
                          mu40dims::bar_lenz / 2 + mu40dims::layer_wallthick / 2), // offset it by corresponding bar width and length
            al_caseLV,                                                             // its logical volume
            "al_case",                                                             // its name
            layer_veto_floorLV,                                                               // its mother  volume
            false,                                                                 // no boolean operation
            0,                                                                     // copy number
            fCheckOverlaps);                                                       // checking overlaps
        auto al_case2PV = new G4PVPlacement(
            0, // no rotation
            G4ThreeVector(0,
                          0,
                          -(mu40dims::bar_lenz / 2 + mu40dims::layer_wallthick / 2)), // offset it by corresponding bar width and length
            al_caseLV,                                                                // its logical volume
            "al_case",                                                                // its name
            layer_veto_floorLV,                                                                  // its mother  volume
            false,                                                                    // no boolean operation
            1,                                                                        // copy number
            fCheckOverlaps);                                                          // checking overlaps            

        G4RotationMatrix transform_rotate_around_z;
        transform_rotate_around_z.rotateZ(90 * deg);
        // First layer
        auto layer_veto_floor1PV = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(),                                 // rotation
                          G4ThreeVector(0, 0, -mu40dims::module_lgap * 0.5)), // offset),
            layer_veto_floorLV,                                               // its logical volume
            "Layer 0",                                                        // its name
            module_veto_floorLV,                                              // its mother volume
            false,                                                            // no boolean operation
            0,                                                                // copy number (detector number within the world)
            fCheckOverlaps);                                                  // checking overlaps
        // Second layer
        auto layer_veto_floor2PV = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(),                                // rotation
                          G4ThreeVector(0, 0, mu40dims::module_lgap * 0.5)), // offset),
            layer_veto_floorLV,                                              // its logical volume
            "Layer 1",                                                       // its name
            module_veto_floorLV,                                             // its mother volume
            false,                                                           // no boolean operation
            1,                                                               // copy number (detector number within the world)
            fCheckOverlaps);                                                 // checking overlaps
    }

    G4LogicalVolume *Mathusla40_Builder::ConstructEnvironment(G4LogicalVolume *_worldLV)
    {
        // Make the detector box again, so that we can subtract it from earth/air
        auto detector = new G4Box("detector", mu40dims::detector_lenx / 2, mu40dims::detector_leny / 2, mu40dims::detector_lenz / 2);

        //
        // Earth - top
        auto earth = new G4Box("earth", mu40dims::world_lenx / 2, mu40dims::world_leny / 2, mu40dims::env_earth_depth_top / 2);
        // Subtract detector from earth (in case we want to excavate a little bit)
        auto earth_excavated =
            new G4SubtractionSolid("earth_excavated",
                                   earth,
                                   detector,
                                   G4Transform3D(G4RotationMatrix(), // rotation
                                                 G4ThreeVector(mu40dims::detector_ground_offset[0],
                                                               mu40dims::detector_ground_offset[1],
                                                               0.5 * mu40dims::env_earth_depth_top + 0.5 * mu40dims::detector_lenz + mu40dims::detector_ground_offset[2])));
        G4LogicalVolume *earthLV = new G4LogicalVolume(
            earth_excavated,     // its solid
            Material::GroundMix, // its material
            "earth");            // its name
        earthLV->SetVisAttributes(Vis::styles["TransparentBrown"]);

        // Place earth in world
        auto earthPV = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(), // rotation
                          G4ThreeVector(0, 0,
                                        -0.5 * mu40dims::env_earth_depth_top)), // offset
            earthLV,                                                            // its logical volume
            "earth",                                                            // its name
            _worldLV,                                                           // its mother volume
            false,                                                              // no boolean operation
            0,                                                                  // copy number (layer number within a module)
            fCheckOverlaps);                                                    // checking overlaps
        (void)earthPV;

        // earth - mid
        auto earth_mid = new G4Box("earth", mu40dims::world_lenx / 2, mu40dims::world_leny / 2, mu40dims::env_earth_depth_mid / 2);
        // Subtract detector from earth (in case we want to excavate a little bit)
        auto earth_excavated_mid =
            new G4SubtractionSolid("earth_excavated",
                                   earth_mid,
                                   detector,
                                   G4Transform3D(G4RotationMatrix(), // rotation
                                                 G4ThreeVector(mu40dims::detector_ground_offset[0],
                                                               mu40dims::detector_ground_offset[1],
                                                               mu40dims::env_earth_depth_top + 0.5 * mu40dims::env_earth_depth_mid + 0.5 * mu40dims::detector_lenz + mu40dims::detector_ground_offset[2])));
        G4LogicalVolume *earthLV_mid = new G4LogicalVolume(
            earth_excavated_mid, // its solid
            Material::GroundMix, // its material
            "earth");            // its name
        earthLV_mid->SetVisAttributes(Vis::styles["TransparentBrown"]);

        // Place earth in world
        auto earthPV_mid = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(), // rotation
                          G4ThreeVector(0, 0,
                                        -mu40dims::env_earth_depth_top - 0.5 * mu40dims::env_earth_depth_mid)), // offset
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
        auto ceiling_out = new G4Box("ceiling_out", mu40dims::env_ceiling_lenx / 2, mu40dims::env_ceiling_leny / 2, mu40dims::env_ceiling_lenz / 2);
        auto ceiling_inside = new G4Box("ceiling_inside", mu40dims::env_ceiling_lenx / 2 - mu40dims::env_ceiling_concrete_thickness,
                                        mu40dims::env_ceiling_leny / 2 - mu40dims::env_ceiling_concrete_thickness,
                                        mu40dims::env_ceiling_lenz / 2 - mu40dims::env_ceiling_concrete_thickness / 2);
        // Subtract detector from earth (in case we want to excavate a little bit)
        auto ceiling =
            new G4SubtractionSolid("ceiling",
                                   ceiling_out,
                                   ceiling_inside,
                                   G4Transform3D(G4RotationMatrix(), // rotation
                                                 G4ThreeVector(0, 0, -1. * mu40dims::env_ceiling_concrete_thickness)));
        G4LogicalVolume *ceilingLV = new G4LogicalVolume(
            ceiling,            // its solid
            Material::Concrete, // its material
            "ceiling");         // its name
        ceilingLV->SetVisAttributes(Vis::styles["TransparentBrown"]);

        // Place ceiling in world
        auto ceilingPV = new G4PVPlacement(
            G4Transform3D(G4RotationMatrix(), // rotation
                          G4ThreeVector(0, 0,
                                        0.5 * mu40dims::env_ceiling_lenz)), // offset
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
