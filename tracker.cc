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
        ("filename", "Digi file to perform reconstruction", cxxopts::value<std::string>())
        ("r,rawfile", "Filename for simulation truth. If provided, it will be merged to the output.", cxxopts::value<std::string>())
        ("c,config", "Filename for configuration.", cxxopts::value<std::string>()->default_value(""))
        ("s,seed", "Seed for random number generator", cxxopts::value<int>()->default_value("-1"))
        ("p,print_progress", "Print progress every `p` events", cxxopts::value<int>()->default_value("1"));
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
    print("   MATHUSLA SIM Tracker, version ", util::VERSION);
    print("      - MATHUSLA Collaboration");
    print("**************************************************************");

    // Configuration

    // Input and output file
    std::filesystem::path input_filename = args["filename"].as<std::string>();
    std::filesystem::path output_filename = input_filename;
    output_filename.replace_filename(output_filename.stem().string() + "_recon" + output_filename.extension().string());
    auto input_reader = std::make_unique<Tracker::TreeReaderDigi>(input_filename.string());
    auto output_writer = std::make_unique<Tracker::TreeWriterRecon>(output_filename.string(), input_reader.get());
    print("Open digi file", input_filename.string());
    print("  Total entries", input_reader->GetEntries());

    // Track and vertex finder
    auto track_finder = Tracker::TrackFinder(false);
    auto vertex_finder = Tracker::VertexFinder(true);

    for (int i = 0; i < input_reader->GetEntries(); i++)
    {
        if ((i + 1) % args["print_progress"].as<int>() == 0)
            print(util::py::f(" Processing {} / {} ", i, input_reader->GetEntries()));

        // Read and process hits
        std::vector<std::unique_ptr<Tracker::DigiHit>> hits_unique = input_reader->GetEntry(i);
        auto hits_groupped = input_reader->ProcessHits(hits_unique);

        // 1. Track finding
        Tracker::TrackList tracks_found;      // Use unique_ptr to auto manage the memory
        std::vector<Tracker::Track *> tracks; // A raw pointer list that duplicates tracks_found
        for (auto &hits_group : hits_groupped)
        {
            std::vector<Tracker::DigiHit *> hits = hits_group.second;
            track_finder.FindAll(hits);
            // auto summary = track_finder.Summary();
            Tracker::TrackList track_found_tmp = track_finder.GetResults();
            for (auto &track : track_found_tmp)
            {
                tracks_found.push_back(std::move(track));
                tracks.push_back(track.get());
            }
        }

        // 2. vertex finding
        vertex_finder.FindAll(tracks);
        Tracker::VertexLilst vertices_found = vertex_finder.GetResults();

        // 3. Append to the output file
        output_writer->ApplyRecon(tracks_found, vertices_found);
        output_writer->Fill();
    }

    output_writer->Write();
    output_writer->Close();

    // Write file to disk
}