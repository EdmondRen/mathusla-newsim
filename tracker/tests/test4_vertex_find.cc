// stdlib
#include <cmath>
#include <chrono>

// ROOT
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "libs/cxxopts.hpp"
#include "test_util.hh"
#include "../vertex_finder.hh"

void setTableauPalette()
{
    // Define RGB values for Tableau colors
    const int nColors = 10;
    double tableauColors[nColors][3] = {
        {31 / 255.0, 119 / 255.0, 180 / 255.0},  // Blue
        {255 / 255.0, 127 / 255.0, 14 / 255.0},  // Orange
        {44 / 255.0, 160 / 255.0, 44 / 255.0},   // Green
        {214 / 255.0, 39 / 255.0, 40 / 255.0},   // Red
        {148 / 255.0, 103 / 255.0, 189 / 255.0}, // Purple
        {140 / 255.0, 86 / 255.0, 75 / 255.0},   // Brown
        {227 / 255.0, 119 / 255.0, 194 / 255.0}, // Pink
        {127 / 255.0, 127 / 255.0, 127 / 255.0}, // Gray
        {188 / 255.0, 189 / 255.0, 34 / 255.0},  // Olive
        {23 / 255.0, 190 / 255.0, 207 / 255.0}   // Cyan
    };

    // Create ROOT colors
    std::vector<int> colorIndices;
    for (int i = 0; i < nColors; ++i)
    {
        int colorIndex = TColor::GetFreeColorIndex();
        new TColor(colorIndex, tableauColors[i][0], tableauColors[i][1], tableauColors[i][2]);
        colorIndices.push_back(colorIndex);
    }

    // Set the palette
    gStyle->SetPalette(colorIndices.size(), &colorIndices[0]);
}

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

    rng.SetSeed(args["seed"].as<int>());

    int N_TRACKS_PER_EVENT = 5;
    // clang-format on

    // Look at one track
    print("Run track/vertex reconstruction on one event.");
    std::vector<Tracker::DigiHit *> hits = genHitsTracks(N_TRACKS_PER_EVENT);
    auto track_finder = Tracker::TrackFinder(hits, false);
    track_finder.FindAll();
    auto summary = track_finder.Summary();
    print(summary);

    // Make a list of track pointers
    Tracker::TrackList track_found = track_finder.GetResults();
    std::vector<Tracker::Track *> tracks;
    for (const auto &track : track_found)
    {
        tracks.push_back(track.get());
    }

    // return 1;
    auto vertex_finder = Tracker::VertexFinder(tracks, true);
    vertex_finder.FindAll();    


    return 0;
}
