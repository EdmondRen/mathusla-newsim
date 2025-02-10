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

int main(int argc, const char *argv[])
{
    // Setup argument format
    // clang-format off
    cxxopts::Options options("./tracker", "Reconstruct the track and vertex in the simulation");
    options.add_options()
        ("h,help", "Print help")
        ("r,rawfile", "Filename for simulation truth. If provided, it will be merged to the output.", cxxopts::value<std::string>()->default_value(""))
        ("R,raw_reduced", "Use this flag to save just the random seed and generator status from the simulation.", cxxopts::value<bool>()->default_value("false"))
        ("c,config", "Filename for configuration.", cxxopts::value<std::string>()->default_value("../macros/tracker/tracker_default.conf"))
        ("k,save_option", "Save only events with track (1) or vertex (2)  or upward vertex (3)", cxxopts::value<int>()->default_value("0"))
        ("s,seed", "Seed for random number generator", cxxopts::value<int>()->default_value("-1"))
        ("p,print_progress", "Print progress every `p` events", cxxopts::value<int>()->default_value("100"))
        ("d,debug_track", "Enable debugging information for track reconstruction", cxxopts::value<bool>()->default_value("false"))
        ("D,debug_vertex", "Enable debugging information for vertex reconstruction", cxxopts::value<bool>()->default_value("false"))
        ("I,event_index", "The index of events to run the tracker, separated by comma. If only two numbers are provide, will treat it as a range.", cxxopts::value<std::vector<int>>()->default_value("0,-1"))
        ("filename", "Digi file to perform reconstruction", cxxopts::value<std::string>());
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
    print("   MATHUSLA SIM reconstruction, version ", util::VERSION);
    print("      - MATHUSLA Collaboration");
    print("      - Authors: ");
    print("**************************************************************");

    // Input and output file
    int save_option = args["save_option"].as<int>();        // 0: all, 1: track:, 2: vertex
    bool save_raw_reduced = args["raw_reduced"].as<bool>(); // save just the seed from simulation file

    std::filesystem::path input_filename = std::filesystem::canonical(args["filename"].as<std::string>());
    std::filesystem::path output_filename = input_filename;
    output_filename.replace_filename(output_filename.stem().string() + "_recon" + output_filename.extension().string());
    print("Processing input file", input_filename.string());
    print("Output will be saved as", output_filename.string());
    auto input_reader = std::make_unique<Tracker::TreeReaderDigi>(input_filename.string());
    auto output_writer = std::make_unique<Tracker::TreeWriterRecon>(output_filename.string(), input_filename.string(), args["rawfile"].as<std::string>(), save_raw_reduced);
    auto start = std::chrono::high_resolution_clock::now();
    auto start_i = std::chrono::high_resolution_clock::now();
    auto stop_i = std::chrono::high_resolution_clock::now();
    print("Open digi file", input_filename.string());
    print("  Total entries", input_reader->GetEntries());

    // Track and vertex finder
    auto track_finder = Tracker::TrackFinder(args["debug_track"].as<bool>());
    auto vertex_finder = Tracker::VertexFinder(args["debug_vertex"].as<bool>());
    auto vertex_track_finder = Tracker::VertexTrackFinder(args["debug_vertex"].as<bool>());
    // Use Configuration file to set up the finder
    auto setup_filename = args["config"].as<std::string>();
    auto parcard = util::io::ParHandler(setup_filename);
    std::map<std::string, double> config = parcard.GetConfig();
    print("  Using config file:", setup_filename);
    track_finder.Config(config);
    vertex_finder.Config(config);
    vertex_track_finder.Config(config);

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

    // Loop selected events
    std::vector<int> inds_to_run;
    auto event_inds = args["event_index"].as<std::vector<int>>();
    if (event_inds.size() == 2)
    {
        int nmin = event_inds[0];
        int nmax = event_inds[1];
        if (nmax == -1)
            nmax = input_reader->GetEntries();
        for (int i = nmin; i < nmax; i++)
            inds_to_run.push_back(i);
    }
    else
        inds_to_run = event_inds;
    info.total_events += inds_to_run.size();

    // for (int i = 0; i < input_reader->GetEntries(); i++)
    for (auto i : inds_to_run)
    {
        if ((i + 1) % args["print_progress"].as<int>() == 0)
        {
            stop_i = std::chrono::high_resolution_clock::now();
            float duration = std::chrono::duration<float>(stop_i - start_i).count();
            print(util::py::f("  Processing {} / {}, {:.2f} sec.", i + 1, inds_to_run.size(), duration));
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
                track->id = tracks_found.size();
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

        // Update the stats
        if (tracks_found.size() > 0)
            info.nevt_with_track += 1;
        if (vertices_found.size() > 0)
            info.nevt_with_vertex += 1;

        // 3. find number of tracklets associated with each vertex
        std::vector<std::unordered_map<int, int>> track_stats_all; // A counter for number of tracks. Pairs of {number_of_hits: number of track}
        int nvertices_upward = 0;
        if (vertices_found.size() > 0)
        {
            if (args["debug_vertex"].as<bool>())
                print("==================Finding tracklets that are consistent with each vertex===================");
            // // Use the vertex with most tracks
            // int ind_maxtrack = -1;
            // int maxtrack = -1;
            // for (auto j = 0; j < static_cast<int>(vertices_found.size()); j++)
            // {
            //     if (static_cast<int>(vertices_found[j]->track_ids.size()) > maxtrack)
            //     {
            //         maxtrack = vertices_found[j]->track_ids.size();
            //         ind_maxtrack = j;
            //     }
            // }

            // Actually, do this on all vertices
            for (auto j = 0; j < static_cast<int>(vertices_found.size()); j++)
            {
                if (!vertices_found[j]->is_downward)
                {
                    nvertices_upward += 1;
                }

                track_stats_all.emplace_back();
                for (auto &hits_group : hits_groupped)
                {

                    std::vector<Tracker::DigiHit *> hits = hits_group.second;
                    vertex_track_finder.FindAll(hits, vertices_found[j].get());
                    vertex_track_finder.Summary();
                    for (auto &pair : vertex_track_finder.track_stats)
                    {
                        track_stats_all.back()[pair.first] += pair.second;
                    }
                }
            }
        }

        // 4. Append to the output file
        // Disable saving if there is no tracks/vertices
        if ((save_option == 1 && tracks_found.size() == 0) ||
            (save_option == 2 && vertices_found.size() == 0) ||
            (save_option == 3 && nvertices_upward == 0))
            continue;
        output_writer->ApplyRecon(tracks_found, vertices_found, track_stats_all, i);
        output_writer->ApplyCopy(i); // Copy digi & raw data
        output_writer->Fill();
    }
    stop_i = std::chrono::high_resolution_clock::now();

    // Write some metadata to the output file

    float duration = std::chrono::duration<float>(stop_i - start).count();
    std::string config_str = parcard.GetString();
    output_writer->meta_ReconstructionConfigStr = config_str;
    output_writer->FillMetadata();

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