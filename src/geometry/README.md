How to get the numbers from the mathusla40 geometry

```c++
float m = 1000;
float cm = 10;
size_t GEO_DEPTH = 4;

// 0: bar. x: along the bar, y: width, z: thickness
// double bar_lenx_real = 2.25 * m;
double bar_lenx_real = 4.5 * m;
double bar_leny_real = 3.5 * cm;
double bar_lenx = bar_lenx_real * round(900 * cm / bar_lenx_real);
double bar_leny = bar_leny_real * round(900 * cm / bar_leny_real / 2) * 2;
double bar_lenz = 1 * cm;
// 1: layer
int layer_Nbars_x = 1;
int layer_Nbars_y = 1;
int layer_Nbars_x_real = round(bar_lenx / bar_lenx_real);
int layer_Nbars_y_real = round(bar_leny / bar_leny_real);
double layer_wallthick = 0.3 * cm;
double layer_hbeam_width = 8 * 2.54 * cm; // Horizontal HSS 8*8*1/2
double layer_hbeam_thick = 0.5 * 2.54 * cm;
double layer_lenx = std::max({bar_lenx * layer_Nbars_x, bar_leny *layer_Nbars_y});
double layer_leny = layer_lenx;
double layer_lenz = layer_hbeam_width * 2 + bar_lenz + layer_wallthick * 2 + 1 * mm;
// 2: tower module
int module_Nlayers = 6;
double module_lgap = 0.8 * m; // Gap between layers
double module_lenx = 10.7 * m;
double module_leny = 10.7 * m;
double module_lenz = module_lgap * (module_Nlayers - 1) + layer_lenz;
std::vector<int> module_layers_zdirection = {kZAxis, kZAxis, kZAxis, kZAxis, kZAxis, kZAxis};
std::vector<int> module_layers_xdirection = {kXAxis, kYAxis, kXAxis, kYAxis, kXAxis, kYAxis};
std::vector<double> module_layers_xoffset = {0, 0, 0, 0, 0}; // From the module x-center
std::vector<double> module_layers_yoffset = {0, 0, 0, 0, 0}; // From the module y-center
std::vector<double> module_layers_zoffset = {layer_lenz / 2,
                                                layer_lenz / 2 + module_lgap * 1,
                                                layer_lenz / 2 + module_lgap * 2,
                                                layer_lenz / 2 + module_lgap * 3,
                                                layer_lenz / 2 + module_lgap * 4,
                                                layer_lenz / 2 + module_lgap * 5}; // From the module z-BOTTOM
double module_vbeam_width = 14 * 2.54 * cm;                                     // Vertical HSS 14*14*5/8
double module_vbeam_thick = 0.625 * 2.54 * cm;
// 3. Entire detector
int detector_Ntowers_x = 4;
int detector_Ntowers_y = 4;
double detector_decay_vol_height = 11 * m;
double detector_lenx = detector_Ntowers_x * module_lenx;
double detector_leny = detector_Ntowers_y * module_leny;
double detector_lenz = module_lenz;
std::vector<double> detector_modules_xoffset = {-1.5 * module_lenx, -0.5 * module_lenx, 0.5 * module_lenx, 1.5 * module_lenx};
std::vector<double> detector_modules_yoffset = detector_modules_xoffset;
std::vector<double> detector_ground_offset = {0, 0, detector_decay_vol_height + module_lgap + layer_lenz}; // x and y offsets are relative to detector box center, z offset is from the bottom.
// Surroundings
double env_earth_depth_top = 0.1 * m;
double env_earth_depth_mid = 100 * m;
double env_air_depth = 40 * m;
double env_ceiling_lenx = detector_lenx + 9.5 * m;
double env_ceiling_leny = detector_lenx + 9.5 * m;
double env_ceiling_lenz = detector_lenz + 1 * m + detector_ground_offset[2]; // about 18 m
double env_ceiling_concrete_thickness = 10 * cm;
double env_floor_iron_thickness = 2 * cm;
// World
double world_lenx = 300 * m;
double world_leny = 200 * m;
double world_lenz = 210 * m;

// Veto layers
// Those are treated differently.
// vf_: veto_floor
double vf_panel_lenx = detector_Ntowers_x * module_lenx;
double vf_panel_leny = detector_Ntowers_y * module_leny;
double vf_zoffset_1 = layer_lenz * 0.5;
double vf_zoffset_2 = layer_lenz * 0.5 + module_lgap;
int vf_Nbars_x_layer1 = ceil(0.5 * vf_panel_lenx / bar_leny_real) * 2;
int vf_Nbars_y_layer1 = ceil(0.5 * vf_panel_leny / bar_lenx_real) * 2;
int vf_Nbars_x_layer2 = ceil(0.5 * vf_panel_lenx / bar_lenx_real) * 2;
int vf_Nbars_y_layer2 = ceil(0.5 * vf_panel_leny / bar_leny_real) * 2;

// vw_: veto_wall
double vw_panel_lenz = detector_decay_vol_height + module_lgap + layer_lenz;
double vw_panel_leny = vf_panel_leny;
int vw_Nbars_z_layer1 = ceil(0.5 * vw_panel_lenz / bar_leny_real) * 2;
int vw_Nbars_y_layer1 = ceil(0.5 * vw_panel_leny / bar_lenx_real) * 2;
int vw_Nbars_z_layer2 = ceil(0.5 * vw_panel_lenz / bar_lenx_real) * 2;
int vw_Nbars_y_layer2 = ceil(0.5 * vw_panel_leny / bar_leny_real) * 2;

cout << "Decay volume bottom offset: "  << detector_ground_offset[2] <<endl;
cout << "Backwall bottom offset: "  << -module_lenx + detector_ground_offset[2] + 75*cm <<endl;

cout << "Veto wall height: "  << vw_panel_lenz <<endl;


```