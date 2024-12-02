// stdlib
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <stdlib.h> // For Ubuntu Linux

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

// Project includes
#include "libs/cxxopts.hpp"
#include "util.hh"
#include "geometry/TestStand_UofT.hh"
#include "geometry/Mathusla40.hh"

// Container for simulation hits
struct SimHit
{
    // Necessary parameters
    int pdg_id; // PDG identifier
    int track_id;
    std::vector<int> copy_numbers;
    float px, py, pz, edep;
    float x, y, z, t;

    // Derived parameters
    //


    u_int64_t det_id;

    SimHit() = default;

    SimHit(float x_position,
           float y_position,
           float z_position,
           float t_position,
           float edep_total,
           float x_momentum,
           float y_momentum,
           float z_momentum,
           int _track_id,
           int _pdg_id,
           std::vector<int> _copy_numbers)
        : x(x_position), y(y_position), z(z_position), t(t_position), edep(edep_total),
          px(x_momentum), py(y_momentum), pz(z_momentum),
          track_id(_track_id), pdg_id(_pdg_id), copy_numbers(_copy_numbers)
    {
        //
        det_id = copy_numbers[0];
        for (size_t i = 1; i < copy_numbers.size(); i++)
        {
            det_id += copy_numbers[i] * std::pow(10, 5 + (i - 1) * 3);
        }
    }
};

// Container for Digi hits
class DigiHit
{
public:
    std::size_t index;
    int direction; // where the bar is pointing to
    double x, ex;
    double y, ey;
    double z, ez;
    double t, et;
    double px, py, pz, edep; // momentum of particle which made the hit
    double particle_mass;
    double particle_energy;
    int pdg_id;
    int track_id;
    int layer_id = 0;
    int type = -1;

    // Internal states
    int xdirection; // xdirection: where the long end is pointing to, one of {0,1,2}
    int zdirection; // zdirection: where the verticle end is pointing to, one of {0,1,2}
    int track_id_min = 9999999999;
    std::vector<SimHit *> hits;

    void AddHit(SimHit *hit, int _direction_)
    {
        hits.push_back(hit);
        if (hit->track_id < track_id_min)
        {
            track_id_min = hit->track_id;
            pdg_id = hit->pdg_id;
            edep = hit->edep;
            px = hit->px;
            py = hit->py;
            pz = hit->pz;
        }
    }
};

struct DigiConfig
{
    float resolution_x, resolution_y, resolution_z, resolution_t;
    float time_limit;
    float sipm_energy_threshold;

    DigiConfig() = default;
    DigiConfig(float _resolution_x, float _resolution_y, float _resolution_z, float _resolution_t, float _time_limit, float _sipm_energy_threshold) : resolution_x(_resolution_x), resolution_y(_resolution_y), resolution_z(_resolution_z), resolution_t(_resolution_t), time_limit(_time_limit), sipm_energy_threshold(_sipm_energy_threshold) {}
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
    std::vector<int> *Hit_copyNumber = nullptr;

    //
    std::vector<std::vector<int>> Hit_copyNumber_split;

    std::string SimulationName;
    std::string Geometry;
    std::string Generator;

    InputTreeHandeler(std::string filename)
    {
        inputFile = TFile::Open(filename.c_str());
        if (!inputFile)
            return;

        treeRaw = (TTree *)inputFile->Get("raw");
        treeMetadata = (TTree *)inputFile->Get("metadata");
        entries = treeRaw->GetEntries();

        // Read metadata
        treeMetadata->SetBranchAddress("SimulationName", &SimulationName);
        treeMetadata->SetBranchAddress("Geometry", &Geometry);
        treeMetadata->SetBranchAddress("Generator", &Generator);
        treeMetadata->GetEntry(0); // Load into buffer

        // Setup tree pointer to data buffer
        treeRaw->SetBranchAddress("Hit_x", &Hit_x);
        treeRaw->SetBranchAddress("Hit_y", &Hit_y);
        treeRaw->SetBranchAddress("Hit_z", &Hit_z);
        treeRaw->SetBranchAddress("Hit_t", &Hit_t);
        treeRaw->SetBranchAddress("Hit_edep", &Hit_edep);
        treeRaw->SetBranchAddress("Hit_px", &Hit_px);
        treeRaw->SetBranchAddress("Hit_py", &Hit_py);
        treeRaw->SetBranchAddress("Hit_pz", &Hit_pz);
        treeRaw->SetBranchAddress("Hit_trackID", &Hit_trackID);
        treeRaw->SetBranchAddress("Hit_pdgID", &Hit_pdgID);
        treeRaw->SetBranchAddress("Hit_copyNumber", &Hit_copyNumber);
        // Speed up reading by turning off all other branches
        treeRaw->SetBranchStatus("Gen_*", 0);
    }

    ~InputTreeHandeler()
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

            //
            int delimiter = -1;
            Hit_copyNumber_split = util::vector::splitVectorByDelimiter(*Hit_copyNumber, delimiter);

            entry_counter += 1;
            for (uint i = 0; i < Hit_x->size(); i++)
            {
                hits.push_back(new SimHit((*Hit_x)[i],
                                          (*Hit_y)[i],
                                          (*Hit_z)[i],
                                          (*Hit_t)[i],
                                          (*Hit_edep)[i],
                                          (*Hit_px)[i],
                                          (*Hit_py)[i],
                                          (*Hit_pz)[i],
                                          (*Hit_trackID)[i],
                                          (*Hit_pdgID)[i],
                                          Hit_copyNumber_split[i]));
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
    std::vector<float> *Digi_x = nullptr;
    std::vector<float> *Digi_y = nullptr;
    std::vector<float> *Digi_z = nullptr;
    std::vector<float> *Digi_t = nullptr;
    std::vector<float> *Digi_edep = nullptr;
    std::vector<float> *Digi_px = nullptr;
    std::vector<float> *Digi_py = nullptr;
    std::vector<float> *Digi_pz = nullptr;
    std::vector<float> *Digi_trackID = nullptr;
    std::vector<double> *Digi_pdgID = nullptr;
    std::vector<int> *Digi_direction = nullptr; // Last three digits of Direction indicates the direction of the bar. For example, 012 means x->x, y->y, z->z
    std::vector<int> *Digi_layerID = nullptr;   // Which layer the hit is from. Layer is obtained from the copy number of depth 1 in GENAT4
    std::vector<int> *Digi_type = nullptr;      // Soure of the event. -1: noise, 0: GUN, 1: PARMA, 2: CRY
    std::vector<int> *Digi_hitInds = nullptr;   // The index of truth hits of each digitized hit

    // Buffer for metadata;
    std::string SimulationName;
    std::string Geometry;
    std::string Generator;
    std::vector<float> Uncertainty;

    OutputTreeHandeler(std::string filename)
    {
        auto output_tree_name = "digi";
        outputFile = TFile::Open(filename.c_str(), "RECREATE");
        outputTreeRaw = new TTree(output_tree_name, "Digitized Tree");
        outputTreeMetadata = new TTree("metadata", "Metadata for digitization");

        // Read metadata
        outputTreeMetadata->Branch("SimulationName", &SimulationName);
        outputTreeMetadata->Branch("Geometry", &Geometry);
        outputTreeMetadata->Branch("Generator", &Generator);
        outputTreeMetadata->Branch("Uncertainty", "std::vector<float>", &Uncertainty);

        // Setup tree pointer to data buffer
        outputTreeRaw->Branch("Digi_x", "std::vector<float>", &Digi_x);
        outputTreeRaw->Branch("Digi_y", "std::vector<float>", &Digi_y);
        outputTreeRaw->Branch("Digi_z", "std::vector<float>", &Digi_z);
        outputTreeRaw->Branch("Digi_t", "std::vector<float>", &Digi_t);
        outputTreeRaw->Branch("Digi_edep", "std::vector<float>", &Digi_edep);
        outputTreeRaw->Branch("Digi_px", "std::vector<float>", &Digi_px);
        outputTreeRaw->Branch("Digi_py", "std::vector<float>", &Digi_py);
        outputTreeRaw->Branch("Digi_pz", "std::vector<float>", &Digi_pz);
        outputTreeRaw->Branch("Digi_trackID", "std::vector<float>", &Digi_trackID);
        outputTreeRaw->Branch("Digi_pdgID", "std::vector<double>", &Digi_pdgID);
        outputTreeRaw->Branch("Digi_layerID", "std::vector<int>", &Digi_layerID);
        outputTreeRaw->Branch("Digi_type", "std::vector<int>", &Digi_type);
        outputTreeRaw->Branch("Digi_hitInds", "std::vector<int>", &Digi_hitInds);
    }

    ~OutputTreeHandeler()
    {
        outputFile->Close();
    }

    void clear()
    {
        Digi_x->clear();
        Digi_y->clear();
        Digi_z->clear();
        Digi_t->clear();
        Digi_edep->clear();
        Digi_px->clear();
        Digi_py->clear();
        Digi_pz->clear();
        Digi_trackID->clear();
        Digi_pdgID->clear();
        Digi_direction->clear();
        Digi_layerID->clear();
        Digi_type->clear();
        Digi_hitInds->clear();
    }

    void Fill()
    {
        outputTreeRaw->Fill();
        clear();
    }
};

bool time_sort(SimHit *hit1, SimHit *hit2) { return (hit1->t < hit2->t); }

std::vector<DigiHit *> Digitize(std::vector<SimHit *> hits, DigiConfig *config)
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
        double x = (current_remaining_hits[0])->x;
        double y = (current_remaining_hits[0])->y;
        double z = (current_remaining_hits[0])->z;

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
                if (hit->t < t0 + config->time_limit)
                {
                    e_sum += hit->edep;
                    used_hits.push_back(hit);
                }
                else
                {
                    unused_hits.push_back(hit);
                }
            }

            if (e_sum > config->sipm_energy_threshold)
            {
                DigiHit *current_digi = new DigiHit();
                // current_digi->det_id = current_id;
                for (auto hit : used_hits)
                {
                    current_digi->AddHit(hit);
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
}

int main(int argc, const char *argv[])
{
    // Setup argument format
    // clang-format off
    cxxopts::Options options("CRY cosmic generator", "CRY cosmic generator with text output");
    options.add_options()
        ("h,help", "Print help")
        ("filename", "ROOT file to digitize", cxxopts::value<std::string>())
        ("n,noise", "Noise rate [avg number per file]. Set to -1 to disable (default).", cxxopts::value<int>()->default_value("-1"))
        ("w,window", "Noise window [ns]", cxxopts::value<float>());
    options.parse_positional({"filename"});
    auto args = options.parse(argc, argv);
    // clang-format on

    // Show help if the user asks for it
    if (args.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0; // Exit the program after showing help
    }

    // Open the input/output file
    std::filesystem::path input_filename = args["filename"].as<std::string>();
    std::filesystem::path output_filename = input_filename;
    output_filename.replace_extension("_digi.root");
    auto infile = new InputTreeHandeler(input_filename.string());
    auto outfile = new InputTreeHandeler(output_filename.string());
}