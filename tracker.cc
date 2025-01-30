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

int main(int argc, const char *argv[])
{
    // Setup argument format
    // clang-format off
    cxxopts::Options options("./tracker", "Reconstruct the track and vertex in the simulation");
    options.add_options()
        ("h,help", "Print help")
        ("r,rawfile", "Filename for simulation truth. If provided, it will be merged to the output.", cxxopts::value<std::string>())
        ("c,config", "Filename for configuration.", cxxopts::value<std::string>()->default_value(""))
        ("s,seed", "Seed for random number generator", cxxopts::value<int>()->default_value("-1"))
        ("p,print_progress", "Print progress every `p` events", cxxopts::value<int>()->default_value("100"))
        ("d,debug_track", "Enable debugging information for track reconstruction", cxxopts::value<bool>()->default_value("false"))
        ("D,debug_vertex", "Enable debugging information for vertex reconstruction", cxxopts::value<bool>()->default_value("false"))
        ("filename", "Digi file to perform reconstruction", cxxopts::value<std::string>());
    options.parse_positional({"filename"});
    options.positional_help("filename"); // Add positional help description
    auto args = options.parse(argc, argv);
    // clang-format on

    // Show help if the user asks for it
    if (args.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0; // Exit the program after showing help
    }

    print("**************************************************************");
    print("   MATHUSLA SIM reconstruction, version ", util::VERSION);
    print("      - MATHUSLA Collaboration");
    print("      - Authors: ");
    print(" ");
    print("**************************************************************");

    // Configuration

    // Input and output file
    std::filesystem::path input_filename = args["filename"].as<std::string>();
    std::filesystem::path output_filename = input_filename;
    output_filename.replace_filename(output_filename.stem().string() + "_recon" + output_filename.extension().string());
    auto input_reader = std::make_unique<Tracker::TreeReaderDigi>(input_filename.string());
    auto output_writer = std::make_unique<Tracker::TreeWriterRecon>(output_filename.string(), input_reader.get());
    auto start = std::chrono::high_resolution_clock::now();
    auto start_i = std::chrono::high_resolution_clock::now();
    auto stop_i = std::chrono::high_resolution_clock::now();
    print("Open digi file", input_filename.string());
    print("  Total entries", input_reader->GetEntries());

    // Track and vertex finder
    auto track_finder = Tracker::TrackFinder(args["debug_track"].as<bool>());
    auto vertex_finder = Tracker::VertexFinder(args["debug_vertex"].as<bool>());
    struct
    {
        int total_events = 0;
        int total_hits = 0;
        int tf_ntracks = 0;
        int tf_nhits = 0;
        int vf_nvertices = 0;
        int vf_ntracks = 0;
        int nevt_with_track = 0;
        int nevt_with_vertex = 0;
    } info;
    info.total_events += input_reader->GetEntries();

    // Loop all events
    for (int i = 0; i < input_reader->GetEntries(); i++)
    {
        if ((i + 1) % args["print_progress"].as<int>() == 0)
        {
            stop_i = std::chrono::high_resolution_clock::now();
            float duration = std::chrono::duration<float>(stop_i - start_i).count();
            print(util::py::f("  Processing {} / {}, {:.2f} sec.", i, input_reader->GetEntries(), duration));
            start_i = stop_i;
        }

        // Read and process hits
        std::vector<std::unique_ptr<Tracker::DigiHit>> hits_unique = input_reader->GetEntry(i);
        std::unordered_map<int, std::vector<Tracker::DigiHit *>> hits_groupped = input_reader->ProcessHits(hits_unique);
        info.total_hits += hits_unique.size();

        // 1. Track finding
        Tracker::TrackList tracks_found, track_found_tmp; // Use unique_ptr to auto manage the memory
        std::vector<Tracker::Track *> tracks;             // A raw pointer list that duplicates tracks_found
        for (auto &hits_group : hits_groupped)
        {
            std::vector<Tracker::DigiHit *> hits = hits_group.second;
            track_finder.FindAll(hits);
            auto summary = track_finder.Summary();
            track_found_tmp = track_finder.GetResults();
            for (auto &track : track_found_tmp)
            {
                tracks_found.push_back(std::move(track));
                tracks.push_back(tracks_found.back().get());
            }

            info.tf_ntracks += track_finder.info_ntracks;
            info.tf_nhits += track_finder.info_nhits;
        }

        // 2. vertex finding
        vertex_finder.FindAll(tracks);
        auto summary = vertex_finder.Summary();
        Tracker::VertexLilst vertices_found = vertex_finder.GetResults();
        info.vf_nvertices += vertex_finder.info_nvertices;
        info.vf_ntracks += vertex_finder.info_ntracks;

        // 3. Append to the output file
        output_writer->ApplyRecon(tracks_found, vertices_found);
        output_writer->Fill();

        // Update the stats
        if (tracks_found.size() > 0)
            info.nevt_with_track += 1;
        if (vertices_found.size() > 0)
            info.nevt_with_vertex += 1;
    }

    float duration = std::chrono::duration<float>(stop_i - start).count();

    print("===============================================");
    print("  Reconstruction Summary:");
    print("    Events/Hits:", info.total_events, "/", info.total_hits);
    print("    Tracks found:", info.tf_ntracks);
    print("       hits used:", info.tf_nhits);
    print("    Vertices found:", info.vf_nvertices);
    print("       tracks used:", info.vf_ntracks);
    print("    Events with track:", info.nevt_with_track);
    print("    Events with vertex:", info.nevt_with_vertex);
    print("    Time used [s]:", duration);
    print("===============================================");

    output_writer->Write();
    output_writer->Close();

    // Write file to disk
}