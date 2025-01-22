// stdlib
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <stdlib.h> // For Ubuntu Linux

// Geant4
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

// Project includes
#include "libs/cxxopts.hpp"
#include "util.hh"
// Include ALL detector geometries that you want to use
#include "geometry/TestStand_UofT.hh"
#include "geometry/Mathusla40.hh"

TRandom3 generator;

std::string util::globals::PROJECT_SOURCE_DIR = "";
std::map<std::string, int> DIGI_TYPE_MAP;

struct DigiConfig
{

    // Speed of light in the fiber
    float c = CLHEP::c_light / 1.89;

    // Coincident timing resolution.
    // == 1/sqrt(2) times single channel resolution
    float time_resolution;

    // Position resolution
    float position_resolution = time_resolution * c;

    // Time and energy cut
    float time_limit;
    float sipm_energy_threshold;

    // Bar dimensions
    float bar_width_x, bar_width_y, bar_width_z;

    DigiConfig() = default;
    DigiConfig(float _time_resolution, float _time_limit, float _sipm_energy_threshold) : time_resolution(_time_resolution), time_limit(_time_limit), sipm_energy_threshold(_sipm_energy_threshold)
    {
        position_resolution = time_resolution * c;
    }
};

// Container for simulation hits
struct SimHit
{
    // Necessary parameters
    int index;
    float x, y, z, t;
    float edep;
    float px, py, pz;
    int track_id;
    int pdg_id; // PDG identifier
    u_int64_t det_id;

    // Derived parameters
    //
    float *pos_vec[3] = {&x, &y, &z};

    SimHit() = default;

    SimHit(int i,
           float x_position,
           float y_position,
           float z_position,
           float t_position,
           float edep_total,
           float x_momentum,
           float y_momentum,
           float z_momentum,
           int _track_id,
           int _pdg_id,
           u_int64_t _det_id)
        : index(i), x(x_position), y(y_position), z(z_position), t(t_position), edep(edep_total),
          px(x_momentum), py(y_momentum), pz(z_momentum),
          track_id(_track_id), pdg_id(_pdg_id), det_id(_det_id) {}
};

// Container for Digi hits
class DigiHit
{
public:
    std::size_t index;
    int direction; // where the bar is pointing to
    double x, y, z, t;
    double px, py, pz, edep; // momentum of particle which made the hit
    double particle_mass;
    double particle_energy;
    int pdg_id;
    int track_id = 999999999;
    unsigned long long detector_id = 0;
    int type = -1;

    // Internal states
    std::vector<SimHit *> hits;

    G4ThreeVector ydirection; // xdirection: where the long end is pointing to, one of {0,1,2}
    G4ThreeVector zdirection; // zdirection: where the verticle end is pointing to, one of {0,1,2}
    long x_dirc_ind = 0;
    long y_dirc_ind = 0;
    long z_dirc_ind = 0;

    double *pos_vec[3] = {&x, &y, &z};

    void AddHit(SimHit *hit, MuGeoBuilder::BarPosition bar_position)
    {
        hits.push_back(hit);
        if (hit->track_id < track_id)
        {
            track_id = hit->track_id;
            pdg_id = hit->pdg_id;
            detector_id = hit->det_id;
            float hit_pos[3] = {hit->x, hit->y, hit->z};

            // Find the y direction index
            for (y_dirc_ind = 0; y_dirc_ind < 3; y_dirc_ind++)
            {
                if (round(abs(bar_position.y_side_direction[y_dirc_ind])) == 1)
                    break;
            }

            // Find the z direction index
            for (z_dirc_ind = 0; z_dirc_ind < 3; z_dirc_ind++)
            {
                if (round(abs(bar_position.z_side_direction[z_dirc_ind])) == 1)
                    break;
            }

            // Find the x direction index, which is the other direction
            x_dirc_ind = 3 - abs(y_dirc_ind) - abs(z_dirc_ind);
            // print("yind, zind", y_dirc_ind, z_dirc_ind);
            // print("y_dir, z_dir", bar_position.y_side_direction.x(), bar_position.y_side_direction.y(), bar_position.y_side_direction.z(), bar_position.z_side_direction.x(), bar_position.z_side_direction.y(), bar_position.z_side_direction.z());

            *pos_vec[y_dirc_ind] = bar_position.bar_center_coord[y_dirc_ind];
            *pos_vec[z_dirc_ind] = bar_position.bar_center_coord[z_dirc_ind];
            *pos_vec[x_dirc_ind] = hit_pos[x_dirc_ind];

            direction = x_dirc_ind * 100 + y_dirc_ind * 10 + z_dirc_ind;
        }
    }

    void Digitize(DigiConfig &config)
    {
        double e_sum = 0;
        double long_direction_sum = 0.0;
        double t_sum = 0;

        for (auto hit : this->hits)
        {
            e_sum += hit->edep;
            t_sum += hit->t * hit->edep;
            long_direction_sum += (*hit->pos_vec[this->x_dirc_ind]) * hit->edep;
        }

        this->edep = e_sum;
        this->t = t_sum / e_sum;
        *this->pos_vec[x_dirc_ind] = long_direction_sum / e_sum;

        // TIME AND POSITION SMEARING
        // we see the random number generator with a number that should be completly random:
        // the clock time times the layer index times the number of digis

        // Time smearing
        this->t += generator.Gaus(0.0, config.time_resolution);
        *this->pos_vec[x_dirc_ind] += generator.Gaus(0.0, config.position_resolution);
    }

    void DigitizeNoise(MuGeoBuilder::BarPosition bar_position, DigiConfig &config, float time_window)
    {

        // Make a fake hit
        auto hit = new SimHit();
        hit->index = 9999999;
        hits.push_back(hit);


        track_id = -1;
        pdg_id = -1;
        detector_id = -1;
        type = -1;

        // Find the y direction index
        for (y_dirc_ind = 0; y_dirc_ind < 3; y_dirc_ind++)
        {
            if (round(abs(bar_position.y_side_direction[y_dirc_ind])) == 1)
                break;
        }

        // Find the z direction index
        for (z_dirc_ind = 0; z_dirc_ind < 3; z_dirc_ind++)
        {
            if (round(abs(bar_position.z_side_direction[z_dirc_ind])) == 1)
                break;
        }

        // Find the x direction index, which is the other direction
        x_dirc_ind = 3 - abs(y_dirc_ind) - abs(z_dirc_ind);

        // Set the results
        *this->pos_vec[y_dirc_ind] = bar_position.bar_center_coord[y_dirc_ind];
        *this->pos_vec[z_dirc_ind] = bar_position.bar_center_coord[z_dirc_ind];
        *this->pos_vec[x_dirc_ind] = bar_position.bar_center_coord[x_dirc_ind] + generator.Uniform(config.bar_width_x) - config.bar_width_x * 0.5;

        this->t = generator.Uniform(time_window * CLHEP::s * 2) - time_window;
        this->direction = x_dirc_ind * 100 + y_dirc_ind * 10 + z_dirc_ind;
    }
};

// --------------------------------------------------------
// Handeler of the simulation data
// * Load event by event into a vector of SimHit
class InputTreeHandeler
{
public:
    TFile *inputFile;
    TTree *treeRaw;
    TTree *treeMetadata;
    int entries;
    int entry_counter = 0;

    std::vector<SimHit *> hits;
    std::vector<float> *Hit_x = nullptr;
    std::vector<float> *Hit_y = nullptr;
    std::vector<float> *Hit_z = nullptr;
    std::vector<float> *Hit_t = nullptr;
    std::vector<float> *Hit_edep = nullptr;
    std::vector<float> *Hit_px = nullptr;
    std::vector<float> *Hit_py = nullptr;
    std::vector<float> *Hit_pz = nullptr;
    std::vector<int> *Hit_trackID = nullptr;
    std::vector<double> *Hit_pdgID = nullptr;
    std::vector<double> *Hit_detectorID = nullptr;

    std::string SimulationName;
    std::string Geometry;
    std::string Generator;
    char SimulationName_buf[200];
    char Geometry_buf[200];
    char Generator_buf[200];

    InputTreeHandeler(std::string filename)
    {
        inputFile = TFile::Open(filename.c_str());
        if (!inputFile)
            return;

        treeRaw = (TTree *)inputFile->Get("raw");
        treeMetadata = (TTree *)inputFile->Get("metadata");
        entries = treeRaw->GetEntries();

        // Read metadata
        treeMetadata->SetBranchAddress("SimulationName", SimulationName_buf);
        treeMetadata->SetBranchAddress("Geometry", Geometry_buf);
        treeMetadata->SetBranchAddress("Generator", Generator_buf);
        treeMetadata->GetEntry(0);           // Load into buffer
        SimulationName = SimulationName_buf; // Turn char[] into string
        Geometry = Geometry_buf;             // Turn char[] into string
        Generator = Generator_buf;           // Turn char[] into string
        SimulationName.erase(SimulationName.find_last_not_of("\n") + 1);
        Geometry.erase(Geometry.find_last_not_of("\n") + 1);
        Generator.erase(Generator.find_last_not_of("\n") + 1);

        // Setup tree pointer to data buffer
        treeRaw->SetBranchAddress("Hit_x", &Hit_x);
        treeRaw->SetBranchAddress("Hit_y", &Hit_y);
        treeRaw->SetBranchAddress("Hit_z", &Hit_z);
        treeRaw->SetBranchAddress("Hit_t", &Hit_t);
        treeRaw->SetBranchAddress("Hit_px", &Hit_px);
        treeRaw->SetBranchAddress("Hit_py", &Hit_py);
        treeRaw->SetBranchAddress("Hit_pz", &Hit_pz);
        treeRaw->SetBranchAddress("Hit_edep", &Hit_edep);
        treeRaw->SetBranchAddress("Hit_trackID", &Hit_trackID);
        treeRaw->SetBranchAddress("Hit_pdgID", &Hit_pdgID);
        treeRaw->SetBranchAddress("Hit_detectorID", &Hit_detectorID);
        // Speed up reading by turning off all other branches
        treeRaw->SetBranchStatus("Gen_*", 0);
    }

    ~InputTreeHandeler()
    {
        inputFile->Close();
    }

    void Close()
    {
        inputFile->Close();
    }

    void Load()
    {
        // Clear previous hits first
        for (auto hit : hits)
            delete hit;
        hits.clear();

        if (entry_counter < entries)
        {
            treeRaw->GetEntry(entry_counter);

            // int delimiter = -1;
            // Hit_detectorID_split = util::vector::splitVectorByDelimiter(*Hit_detectorID, delimiter);

            entry_counter += 1;
            for (uint i = 0; i < Hit_x->size(); i++)
            {
                hits.push_back(new SimHit(i,
                                          (*Hit_x)[i],
                                          (*Hit_y)[i],
                                          (*Hit_z)[i],
                                          (*Hit_t)[i],
                                          (*Hit_edep)[i],
                                          (*Hit_px)[i],
                                          (*Hit_py)[i],
                                          (*Hit_pz)[i],
                                          (*Hit_trackID)[i],
                                          (*Hit_pdgID)[i],
                                          (*Hit_detectorID)[i]));
            }
        }
    }
};

// --------------------------------------------------
// Handeler of the digitized data
class OutputTreeHandeler
{
public:
    TFile *outputFile;
    TTree *outputTreeRaw;
    TTree *outputTreeMetadata;

    // Buffer for output data
    std::vector<float> Digi_x;
    std::vector<float> Digi_y;
    std::vector<float> Digi_z;
    std::vector<float> Digi_t;
    std::vector<float> Digi_edep;
    std::vector<int> Digi_trackID;
    std::vector<int> Digi_pdgID;
    std::vector<long long> Digi_detectorID; // Which layer the hit is from. Layer is obtained from the copy number of depth 1 in GENAT4
    std::vector<int> Digi_type;             // Soure of the event. -1: noise, 0: GUN, 1: PARMA, 2: CRY
    std::vector<int> Digi_hitInds;          // The index of truth hits of each digitized hit
    std::vector<int> Digi_direction;        // Indicates the direction of the bar with last three digits . For example, .....012 means x->x, y->y, z->z

    // Buffer for metadata;
    std::string SimulationName;
    std::string Geometry;
    std::string Generator;
    float Uncertainty_t;
    float Uncertainty_x;
    float Uncertainty_y;
    float Uncertainty_z;

    OutputTreeHandeler(std::string filename)
    {
        auto output_tree_name = "digi";
        outputFile = TFile::Open(filename.c_str(), "RECREATE");
        outputTreeRaw = new TTree(output_tree_name, "Digitized Tree");
        outputTreeMetadata = new TTree("metadata", "Metadata for digitization");

        // Write metadata
        outputTreeMetadata->Branch("SimulationName", &SimulationName);
        outputTreeMetadata->Branch("Geometry", &Geometry);
        outputTreeMetadata->Branch("Generator", &Generator);
        outputTreeMetadata->Branch("Uncertainty_t", &Uncertainty_t);
        outputTreeMetadata->Branch("Uncertainty_x", &Uncertainty_x);
        outputTreeMetadata->Branch("Uncertainty_y", &Uncertainty_y);
        outputTreeMetadata->Branch("Uncertainty_z", &Uncertainty_z);

        // Setup tree pointer to data buffer
        outputTreeRaw->Branch("Digi_x", &Digi_x);
        outputTreeRaw->Branch("Digi_y", &Digi_y);
        outputTreeRaw->Branch("Digi_z", &Digi_z);
        outputTreeRaw->Branch("Digi_t", &Digi_t);
        outputTreeRaw->Branch("Digi_edep", &Digi_edep);
        outputTreeRaw->Branch("Digi_trackID", &Digi_trackID);
        outputTreeRaw->Branch("Digi_pdgID", &Digi_pdgID);
        outputTreeRaw->Branch("Digi_detectorID", &Digi_detectorID);
        outputTreeRaw->Branch("Digi_type", &Digi_type);
        outputTreeRaw->Branch("Digi_hitInds", &Digi_hitInds);
        outputTreeRaw->Branch("Digi_direction", &Digi_direction);
    }

    ~OutputTreeHandeler()
    {
        outputFile->Close();
    }

    void clear()
    {
        Digi_x.clear();
        Digi_y.clear();
        Digi_z.clear();
        Digi_t.clear();
        Digi_edep.clear();
        Digi_trackID.clear();
        Digi_pdgID.clear();
        Digi_direction.clear();
        Digi_detectorID.clear();
        Digi_type.clear();
        Digi_hitInds.clear();
    }

    void ExportDigis(std::vector<DigiHit *> digi_list)
    {
        this->clear();

        for (auto digi : digi_list)
        {
            Digi_t.push_back(digi->t);
            Digi_x.push_back(digi->x);
            Digi_y.push_back(digi->y);
            Digi_z.push_back(digi->z);
            Digi_edep.push_back(digi->edep);
            Digi_trackID.push_back(digi->track_id);
            Digi_pdgID.push_back(digi->pdg_id);
            Digi_direction.push_back(digi->direction);
            Digi_detectorID.push_back(digi->detector_id);
            Digi_type.push_back(digi->type);
            for (auto hit : digi->hits)
            {
                Digi_hitInds.push_back(hit->index);
            }
            Digi_hitInds.push_back(-1);
        }
    }

    void Fill()
    {
        outputTreeRaw->Fill();
        clear();
    }

    void Write()
    {
        outputTreeRaw->Write();
    }

    void WriteConfig(DigiConfig &config)
    {
        Uncertainty_t = config.time_resolution;
        Uncertainty_x = config.position_resolution;
        Uncertainty_y = config.bar_width_y / std::sqrt(12);
        Uncertainty_z = config.bar_width_z / std::sqrt(12);
        outputTreeMetadata->Fill();
        outputTreeMetadata->Write();
    }

    void Close()
    {
        outputFile->Close();
    }
};

bool time_sort(SimHit *hit1, SimHit *hit2) { return (hit1->t < hit2->t); }

std::vector<DigiHit *> Digitize(std::vector<SimHit *> hits, DigiConfig &config, MuGeoBuilder::Builder *geobulder)
{
    // this is the vector of digi_hits we will return at the end of the function
    std::vector<DigiHit *> digis;
    // looping through each detector ID
    std::vector<SimHit *> current_hits;
    std::vector<SimHit *> current_remaining_hits = hits;
    std::vector<SimHit *> next_remaining_hits;

    while (current_remaining_hits.size() > 0)
    {
        // current detector id which we are working in
        auto current_id = (current_remaining_hits[0])->det_id;
        // double x = (current_remaining_hits[0])->x;

        // taking out all hits with the same detector id to be digitized, leaving the remaing for the next iteration
        for (auto hit : current_remaining_hits)
        {
            if (hit->det_id == current_id)
            {
                current_hits.push_back(hit);
            }
            else
            {
                next_remaining_hits.push_back(hit);
            }
        }

        // time sorting current hits
        std::sort(current_hits.begin(), current_hits.end(), &time_sort);

        // going through all hits until they are either all added to digis, or dropped
        while (current_hits.size() > 0)
        {

            std::vector<SimHit *> used_hits;
            std::vector<SimHit *> unused_hits;

            double t0 = (current_hits[0])->t;
            double e_sum = 0;

            for (auto hit : current_hits)
            {
                if (hit->t < t0 + config.time_limit)
                {
                    e_sum += hit->edep;
                    used_hits.push_back(hit);
                }
                else
                {
                    unused_hits.push_back(hit);
                }
            }

            // print("esum", e_sum, "len used hit", used_hits.size());

            if (e_sum > config.sipm_energy_threshold)
            {
                DigiHit *current_digi = new DigiHit();
                for (auto hit : used_hits)
                {
                    current_digi->AddHit(hit, geobulder->GetBarPosition(hit->det_id));
                }
                current_digi->index = (digis.size());
                digis.push_back(current_digi);
                current_hits = unused_hits;
            }
            else
            {
                current_hits.erase(current_hits.begin());
            }

        } // while (current_hits.size() > 0)

        // resetting all the sorting vectors, and assigning the next remianing hits to the next iteration for current remaining
        current_remaining_hits.clear();
        current_remaining_hits = next_remaining_hits;
        next_remaining_hits.clear();
        current_hits.clear();
    } // while (current_remaining_hits.size() > 0)

    // At this point, all of the digi_hits in the digi_vector have the hits which will make them up. However, they don't have any of their energy, position, or timing information added.
    // Below, we compute the energy, time, and position of all of the digi hits
    // We incoorporate the time and position smearing into this calculation as well

    // setting digi indices
    int k = 0;
    for (auto digi : digis)
    {
        digi->Digitize(config);
        digi->index = k++;
    }

    return digis;
}

class NoiseMaker
{
public:
    MuGeoBuilder::BarPositionMap bar_map;
    std::vector<unsigned long long> bar_index;
    float noise_window, noise_counts_mean;
    DigiConfig &config;

    // noise_rate: [Hz], noise rate of each bar
    // noise_window: [s], noise will be sampled between [-noise_window, noise_window]
    NoiseMaker(float noise_rate, float _noise_window, DigiConfig &_config, MuGeoBuilder::Builder *geobuilder) : noise_window(_noise_window), config(_config)
    {
        // Get a list of all bars
        this->bar_map = geobuilder->GetBarPositionMap();

        // Give each bar an index starting from 0
        for (const auto &[key, value] : this->bar_map)
        {
            bar_index.push_back(key);
        }

        // Average total number of noise events in all bars
        this->noise_counts_mean = noise_rate * noise_window * bar_index.size() * 2;
    }

    std::vector<DigiHit *> run()
    {
        std::vector<DigiHit *> digis; // this is the vector of digi_hits we will return at the end of the function
        std::vector<int> bars_selected;

        if (noise_counts_mean <= 0)
            return digis;

        int noise_counts_this = generator.Poisson(noise_counts_mean);
        for (int i = 0; i < noise_counts_this; i++)
        {
            int bar_ind = generator.Integer(bar_index.size());
            auto detector_id = bar_index.at(bar_ind);
            auto bar_position = bar_map.at(detector_id);

            DigiHit *current_digi = new DigiHit();
            current_digi->DigitizeNoise(bar_position, config, noise_window);
            current_digi->index = i;
            current_digi->detector_id = detector_id;
            digis.push_back(current_digi);
        }
        return digis;
    }
};

int main(int argc, const char *argv[])
{
    // Setup argument format
    // clang-format off
    cxxopts::Options options("./digitizer", "Digitize the MATHUSLA simulation events");
    options.add_options()
        ("h,help", "Print help")
        ("filename", "ROOT file to digitize", cxxopts::value<std::string>())
        ("s,seed", "Seed for random number generator", cxxopts::value<int>()->default_value("-1"))
        ("t,time_resolution", "Coincidence time resolution [ns].", cxxopts::value<float>()->default_value("1"))
        ("T,time_limit", "Time limit [ns]", cxxopts::value<float>()->default_value("20"))
        ("E,energy_threshold", "Energy threshold for a digi [MeV]", cxxopts::value<float>()->default_value("0.65"))
        ("p,print_progress", "Print progress every `p` events", cxxopts::value<int>()->default_value("1"))
        ("n,noise_rate", "Noise rate [avg number per file]. Set to -1 to disable (default).", cxxopts::value<float>()->default_value("-1"))
        ("w,noise_window", "Noise window [s], noise will be sampled between [-noise_window, noise_window]", cxxopts::value<float>()->default_value("1000e-9"));
    options.parse_positional({"filename"});
    auto args = options.parse(argc, argv);
    // clang-format on

    // Show help if the user asks for it
    if (args.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0; // Exit the program after showing help
    }

    print("**************************************************************");
    print("   MATHUSLA SIM Digitizer, version ", util::VERSION);
    print("      - MATHUSLA Collaboration");
    print("**************************************************************");

    // Make a map between generator name and the Digi_type
    DIGI_TYPE_MAP.insert({"gun", 0});
    DIGI_TYPE_MAP.insert({"cry", 1});
    DIGI_TYPE_MAP.insert({"parma", 2});
    DIGI_TYPE_MAP.insert({"filereader", 3});
    DIGI_TYPE_MAP.insert({"noise", -1});

    // Setup random number generator
    generator.SetSeed(args["seed"].as<int>());
    auto print_progress = args["print_progress"].as<int>();

    // make a configuration for digitizer
    float _time_resolution = args["time_resolution"].as<float>();        // [ns]
    float _time_limit = args["time_limit"].as<float>();                  // [ns]
    float _sipm_energy_threshold = args["energy_threshold"].as<float>(); // [MeV]
    auto config = DigiConfig(_time_resolution, _time_limit, _sipm_energy_threshold);

    // Open the input/output file
    std::filesystem::path input_filename = args["filename"].as<std::string>();
    std::filesystem::path output_filename = input_filename;
    output_filename.replace_extension("digi.root");
    auto infile = new InputTreeHandeler(input_filename.string());
    auto outfile = new OutputTreeHandeler(output_filename.string());

    // Build all geometries and select the one based on metadata of infile
    std::unordered_map<std::string, MuGeoBuilder::Builder *> _det_map_;
    _det_map_["uoft1"] = new MuGeoBuilder::Uoft1_Builder();
    _det_map_["mu40v0"] = new MuGeoBuilder::Mathusla40_Builder();
    auto _det_selected_ = _det_map_[infile->Geometry];
    _det_selected_->Construct();

    print("Digitizer > Building Select geometry: ", infile->Geometry, "**");
    print("Digitizer > Finished building geometry ");
    print("Digitizer > Running");

    int digi_type = DIGI_TYPE_MAP[infile->Generator];

    // Write the metatdata into output file
    // 1. Copy the information from input file
    outfile->SimulationName = infile->SimulationName;
    outfile->Geometry = infile->Geometry;
    outfile->Generator = infile->Generator;
    // Read geometry information from the GeoBuilder
    if (infile->Geometry == "uoft1")
    {
        config.bar_width_x = MuGeoBuilder::uoftdims::bar_lenx;
        config.bar_width_y = MuGeoBuilder::uoftdims::bar_leny;
        config.bar_width_z = MuGeoBuilder::uoftdims::bar_lenz;
    }
    if (infile->Geometry == "mu40v0")
    {
        config.bar_width_x = MuGeoBuilder::mu40dims::bar_lenx_real;
        config.bar_width_y = MuGeoBuilder::mu40dims::bar_leny_real;
        config.bar_width_z = MuGeoBuilder::mu40dims::bar_lenz;
    }    
    // 2. Add info from config
    outfile->WriteConfig(config);

    // Start noise maker
    auto noise_maker = NoiseMaker(args["noise_rate"].as<float>(), args["noise_window"].as<float>(), config, _det_selected_);

    for (int entry = 0; entry < (infile->entries); entry++)
    {
        // Print progress
        if (entry % print_progress == 0)
            print("Digitizer > ---> Begin of event:", entry);

        // Load hits
        infile->Load();

        // Digitize
        auto digis = Digitize(infile->hits, config, _det_selected_);
        for (auto digi : digis)
        {
            digi->type = digi_type;
            // Debug, Inspect the digits
            // print("Digi x", digi->x);
        }

        // Make noise hits
        if (args["noise_rate"].as<float>() > 0)
        {
            auto noise_digis = noise_maker.run();
            int last_digi_index = digis.size() > 0 ? digis.back()->index : -1;

            // Set the index and type of the noise hits, then add to list
            for (auto digi : noise_digis)
            {
                digi->type = DIGI_TYPE_MAP["noise"];
                digi->index = digi->index + last_digi_index + 1;
                digis.push_back(digi);
            }
        }

        // Export digis to output file
        outfile->ExportDigis(digis);
        outfile->Fill();

        // Clear the digis
        for (auto digi : digis)
            delete digi;
        

    }

    outfile->Write();
    outfile->Close();
    infile->Close();

    print("Digitizer > Finished. File saved as:", output_filename);
}