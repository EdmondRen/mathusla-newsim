// stdlib
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <chrono>
#include <stdlib.h> // For Ubuntu Linux

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

// Project includes
#include "libs/cxxopts.hpp"
#include "util.hh"
#include "types_and_helper.hh"
#include "track_finder.hh"
#include "vertex_finder.hh"
// Global variables
#include "util_globals.hh"

std::string util::globals::PROJECT_SOURCE_DIR = "";

class ReconFileHandeler
{
public:
    TFile *outputFile;
    TTree *outputTree;
    TFile *reconFile;
    TTree *reconTree;
    Long64_t entry_current;

    std::vector<std::unique_ptr<Tracker::Track>> tracks;
    std::vector<std::unique_ptr<Tracker::Vertex>> vertices;

    // Pointers to branches
    std::vector<float> *Track_x0;
    std::vector<float> *Track_y0;
    std::vector<float> *Track_z0;
    std::vector<float> *Track_t0;
    std::vector<float> *Track_kx;
    std::vector<float> *Track_ky;
    std::vector<float> *Track_kz;
    std::vector<float> *Track_kt;
    std::vector<float> *Track_cov;
    std::vector<float> *Track_chi2;
    std::vector<int> *Track_id;
    std::vector<int> *Track_iv_ind;
    std::vector<int> *Track_iv_err;
    std::vector<int> *Track_digiInds;
    std::vector<float> *Vertex_x0;
    std::vector<float> *Vertex_y0;
    std::vector<float> *Vertex_z0;
    std::vector<float> *Vertex_t0;
    std::vector<float> *Vertex_cov;
    std::vector<float> *Vertex_chi2;
    std::vector<int> *Vertex_id;
    std::vector<int> *Vertex_trackInds;
    std::vector<int> *Vertex_tracklet_n0;
    std::vector<int> *Vertex_tracklet_n2;
    std::vector<int> *Vertex_tracklet_n3;
    std::vector<int> *Vertex_tracklet_n4p;

    ReconFileHandeler(std::string filename_recon,
                      std::string filename_out)
    {
        auto tree_name = "data";
        reconFile = TFile::Open(filename_recon.c_str());
        reconTree = (TTree *)reconFile->Get(tree_name);

        outputFile = TFile::Open(filename_out.c_str(), "RECREATE");
        outputTree = new TTree(tree_name, "Reconstruction Tree Skimmed");

        // Setup for the address of all branches
        // They could then be accessed by the key names
        // this->Setup(reconTree, outputTree);

        // Track_x0 = this->GetFloatV("Track_x0");
        // Track_y0 = this->GetFloatV("Track_y0");
        // Track_z0 = this->GetFloatV("Track_z0");
        // Track_t0 = this->GetFloatV("Track_t0");
        // Track_kx = this->GetFloatV("Track_kx");
        // Track_ky = this->GetFloatV("Track_ky");
        // Track_kz = this->GetFloatV("Track_kz");
        // Track_kt = this->GetFloatV("Track_kt");
        // Track_cov = this->GetFloatV("Track_cov");
        // Track_chi2 = this->GetFloatV("Track_chi2");
        // Track_id = this->GetIntV("Track_id");
        // Track_iv_ind = this->GetIntV("Track_iv_ind");
        // Track_iv_err = this->GetIntV("Track_iv_err");
        // Track_digiInds = this->GetIntV("Track_digiInds");
        // Vertex_x0 = this->GetFloatV("Vertex_x0");
        // Vertex_y0 = this->GetFloatV("Vertex_y0");
        // Vertex_z0 = this->GetFloatV("Vertex_z0");
        // Vertex_t0 = this->GetFloatV("Vertex_t0");
        // Vertex_cov = this->GetFloatV("Vertex_cov");
        // Vertex_chi2 = this->GetFloatV("Vertex_chi2");
        // Vertex_id = this->GetIntV("Vertex_id");
        // Vertex_trackInds = this->GetIntV("Vertex_trackInds");
        // Vertex_tracklet_n0 = this->GetIntV("Vertex_tracklet_n0");
        // Vertex_tracklet_n2 = this->GetIntV("Vertex_tracklet_n2");
        // Vertex_tracklet_n3 = this->GetIntV("Vertex_tracklet_n3");
        // Vertex_tracklet_n4p = this->GetIntV("Vertex_tracklet_n4p");

        // Setup tree pointer to data buffer
        reconTree->SetBranchAddress("Track_x0", &Track_x0);
        reconTree->SetBranchAddress("Track_y0", &Track_y0);
        reconTree->SetBranchAddress("Track_z0", &Track_z0);
        reconTree->SetBranchAddress("Track_t0", &Track_t0);
        reconTree->SetBranchAddress("Track_kx", &Track_kx);
        reconTree->SetBranchAddress("Track_ky", &Track_ky);
        reconTree->SetBranchAddress("Track_kz", &Track_kz);
        reconTree->SetBranchAddress("Track_kt", &Track_kt);
        reconTree->SetBranchAddress("Track_cov", &Track_cov); // Have to be flattened, each track takes 6x6=36 elements
        reconTree->SetBranchAddress("Track_chi2", &Track_chi2);
        reconTree->SetBranchAddress("Track_id", &Track_id);
        reconTree->SetBranchAddress("Track_iv_ind", &Track_iv_ind);
        reconTree->SetBranchAddress("Track_iv_err", &Track_iv_err);
        reconTree->SetBranchAddress("Track_digiInds", &Track_digiInds);
        reconTree->SetBranchAddress("Vertex_x0", &Vertex_x0);
        reconTree->SetBranchAddress("Vertex_y0", &Vertex_y0);
        reconTree->SetBranchAddress("Vertex_z0", &Vertex_z0);
        reconTree->SetBranchAddress("Vertex_t0", &Vertex_t0);
        reconTree->SetBranchAddress("Vertex_cov", &Vertex_cov); // Have to be flattened, each track takes 6x6=36 elements
        reconTree->SetBranchAddress("Vertex_chi2", &Vertex_chi2);
        reconTree->SetBranchAddress("Vertex_id", &Vertex_id);
        reconTree->SetBranchAddress("Vertex_trackInds", &Vertex_trackInds);
        reconTree->SetBranchAddress("Vertex_tracklet_n0", &Vertex_tracklet_n0);
        reconTree->SetBranchAddress("Vertex_tracklet_n2", &Vertex_tracklet_n2);
        reconTree->SetBranchAddress("Vertex_tracklet_n3", &Vertex_tracklet_n3);
        reconTree->SetBranchAddress("Vertex_tracklet_n4p", &Vertex_tracklet_n4p);

        outputTree->Branch("Track_x0", &(*Track_x0));
        outputTree->Branch("Track_y0", &(*Track_y0));
        outputTree->Branch("Track_z0", &(*Track_z0));
        outputTree->Branch("Track_t0", &(*Track_t0));
        outputTree->Branch("Track_kx", &(*Track_kx));
        outputTree->Branch("Track_ky", &(*Track_ky));
        outputTree->Branch("Track_kz", &(*Track_kz));
        outputTree->Branch("Track_kt", &(*Track_kt));
        outputTree->Branch("Track_cov", &(*Track_cov)); // Have to be flattened, each track takes 6x6=36 elements
        outputTree->Branch("Track_chi2", &(*Track_chi2));
        outputTree->Branch("Track_id", &(*Track_id));
        outputTree->Branch("Track_iv_ind", &(*Track_iv_ind));
        outputTree->Branch("Track_iv_err", &(*Track_iv_err));
        outputTree->Branch("Track_digiInds", &(*Track_digiInds));
        outputTree->Branch("Vertex_x0", &(*Vertex_x0));
        outputTree->Branch("Vertex_y0", &(*Vertex_y0));
        outputTree->Branch("Vertex_z0", &(*Vertex_z0));
        outputTree->Branch("Vertex_t0", &(*Vertex_t0));
        outputTree->Branch("Vertex_cov", &(*Vertex_cov)); // Have to be flattened, each track takes 6x6=36 elements
        outputTree->Branch("Vertex_chi2", &(*Vertex_chi2));
        outputTree->Branch("Vertex_id", &(*Vertex_id));
        outputTree->Branch("Vertex_trackInds", &(*Vertex_trackInds));
        outputTree->Branch("Vertex_tracklet_n0", &(*Vertex_tracklet_n0));
        outputTree->Branch("Vertex_tracklet_n2", &(*Vertex_tracklet_n2));
        outputTree->Branch("Vertex_tracklet_n3", &(*Vertex_tracklet_n3));
        outputTree->Branch("Vertex_tracklet_n4p", &(*Vertex_tracklet_n4p));
    }

    ~ReconFileHandeler()
    {
        reconFile->Close();
        outputFile->Close();
    }

    Long64_t GetEntries() { return reconTree->GetEntries(); }
    void Fill() { outputTree->Fill(); }
    void Write() { outputFile->Write(); }
    void Close() { outputFile->Close(); }

    std::vector<std::vector<int>> splitVector(const std::vector<int> &input, int splitValue)
    {
        std::vector<std::vector<int>> result;
        std::vector<int> current;

        for (int value : input)
        {
            if (value == splitValue)
            {
                if (!current.empty())
                {
                    result.push_back(current);
                    current.clear();
                }
            }
            else
            {
                current.push_back(value);
            }
        }

        // Add the last segment if it's not empty
        if (!current.empty())
        {
            result.push_back(current);
        }

        return result;
    }

    // Read an event from reconstruction file
    // Organize them into vertex and track objects
    // Tracks and Vertices will be accessible via "tracks" and "vertices"
    int GetEntry(Long64_t entry)
    {
        int status = reconTree->GetEntry(entry);
        entry_current = entry;

        this->tracks.clear();
        this->vertices.clear();

        // Make tracks
        for (int i = 1; i < (int)Track_x0->size(); i++)
        {
            auto track = std::make_unique<Tracker::Track>();
            track->id = (*Track_id)[i];
            track->iv_index = (*Track_iv_ind)[i];
            track->iv_error = (*Track_iv_err)[i];
            track->params_full = VectorXd(8);
            track->params_full << (*Track_x0)[i], (*Track_y0)[i], (*Track_z0)[i], (*Track_t0)[i], (*Track_kx)[i], (*Track_ky)[i], (*Track_kz)[i], (*Track_kt)[i];
            auto params = Tracker::Helper::removeElement(track->params_full, track->iv_index + 4);
            params = Tracker::Helper::removeElement(params, track->iv_index);
            track->params = params;
            this->tracks.push_back(std::move(track));
        }

        // Make vertices
        std::vector<std::vector<int>> vertex_trackid_split = splitVector((*Vertex_trackInds), -1);
        for (int i = 1; i < (int)Vertex_x0->size(); i++)
        {
            auto vertex = std::make_unique<Tracker::Vertex>();
            vertex->params = Vector4d((*Vertex_x0)[i], (*Vertex_y0)[i], (*Vertex_z0)[i], (*Vertex_t0)[i]);

            // std::vector<float> cov_vec(Vertex_cov.begin() + i * 16, Vertex_cov.begin() + (i + 1) * 16);
            // vertex->cov = MatrixXd(4,4);
            vertex->track_ids = vertex_trackid_split[i];
            this->vertices.push_back(std::move(vertex));
        }

        return status;
    }
};

int main(int argc, const char *argv[])
{
    // Setup argument format
    // clang-format off
    cxxopts::Options options("./skim_vertex", "Select only events that have signal-like vertex");
    options.add_options()
        ("h,help", "Print help")
        ("o,output", "Output filename. If not given, will be automatically set as {input_filename}_skim.root", cxxopts::value<std::string>())
        ("k,save_option", "Save only events with track (1) or vertex (2)  or upward vertex (3). Default is saving all events.", cxxopts::value<int>()->default_value("0"))
        ("p,print_progress", "Print progress every `p` events", cxxopts::value<int>()->default_value("100"))
        ("filename", "Recon. file to work on", cxxopts::value<std::string>());
    options.parse_positional({"filename"});
    options.positional_help("digit_filename"); // Add positional help description
    auto args = options.parse(argc, argv);
    // clang-format on

    // Show help if the user asks for it
    if (args.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0; // Exit the program after showing help
    }

    print("\n**************************************************************");
    print("   MATHUSLA SIM data skimming tool, version ", util::VERSION);
    print("      - MATHUSLA Collaboration");
    print("      - Authors: ");
    print("**************************************************************");

    // Input and output file
    std::filesystem::path input_filename = std::filesystem::canonical(args["filename"].as<std::string>());
    std::filesystem::path output_filename = input_filename;
    if (args.count("output") == 0)
        output_filename.replace_filename(output_filename.stem().string() + "_skim" + output_filename.extension().string());
    else
        output_filename = args["output"].as<std::string>();
    print("Processing input file", input_filename.string());
    print("Output will be saved as", output_filename.string());

    // The data object will make it easier to read from ROOT tree.
    auto data = std::make_unique<ReconFileHandeler>(input_filename.string(), output_filename.string());

    auto start = std::chrono::high_resolution_clock::now();
    auto start_i = std::chrono::high_resolution_clock::now();
    auto stop_i = std::chrono::high_resolution_clock::now();
    print("Open recon file", input_filename.string());
    print("  Total entries", data->GetEntries());

    // Define the decay volume
    float box_halflen_mm = 10700 * 2;
    float box_fiducial_shrink = 300; // Shrink inside by 30cm
    float box_x[2] = {-box_halflen_mm + box_fiducial_shrink, box_halflen_mm - box_fiducial_shrink};
    float box_y[2] = {-box_halflen_mm + box_fiducial_shrink, box_halflen_mm - box_fiducial_shrink};
    float box_z[2] = {1020 + box_fiducial_shrink, 12000 - box_fiducial_shrink};

    // Loop all events
    auto entries = data->GetEntries();
    for (int i = 0; i < entries; i++)
    {
        if ((i + 1) % args["print_progress"].as<int>() == 0)
        {
            stop_i = std::chrono::high_resolution_clock::now();
            float duration = std::chrono::duration<float>(stop_i - start_i).count();
            print(util::py::f("  Processing {} / {}, {:.2f} sec.", i + 1, entries, duration));
            start_i = stop_i;
        }

        // Read i_th event. We should have data.tracks and data.vertices
        data->GetEntry(i);

        // Check all vertices.
        // It takes only one vertex that passed all cuts to keep this event.
        bool passed = false;
        for (auto &v : data->vertices)
        {

            // 0. Vertex is inside the box
            bool is_inside = (v->params(0) > box_x[0] && v->params(0) < box_x[1] &&
                              v->params(1) > box_y[0] && v->params(1) < box_y[1] &&
                              v->params(2) > box_z[0] && v->params(2) < box_z[1]);

            // 1. No inward tracks
            bool is_outward = true;
            for (auto j : v->track_ids)
            {
                auto &t = data->tracks[j];
                bool track_is_outward = (t->iv_index == 2 && t->params_full(7) > 0) ||
                                        (t->iv_index == 0 && t->params_full(7) > 0);
                if (!track_is_outward)
                {
                    is_outward = false;
                    break;
                }
            }
            print(is_inside, is_outward);

            if (passed = (is_inside && is_outward))
                break;
        }

        // Save the event to the output file if passed cuts.
        if (passed)
        {
            data->Fill();
            print("Vertex passed");
        }
    }
}