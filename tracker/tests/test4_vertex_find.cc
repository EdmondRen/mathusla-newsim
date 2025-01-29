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

std::string util::globals::PROJECT_SOURCE_DIR = "";

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
        ("n,ntracks", "N.", cxxopts::value<int>()->default_value("5"))
        ("e,nevents", "N.", cxxopts::value<int>()->default_value("1000"))
        ("w,noise_window", "Noise window [s], noise will be sampled between [-noise_window, noise_window]", cxxopts::value<float>()->default_value("1000e-9"));
    options.parse_positional({"filename"});
    auto args = options.parse(argc, argv);

    rng.SetSeed(args["seed"].as<int>());

    int N_TRACKS_PER_EVENT =args["ntracks"].as<int>();
    // clang-format on

    // Look at one Event
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

    // Run vertex finder
    auto vertex_finder = Tracker::VertexFinder(tracks, true);
    vertex_finder.FindAll();
    // Get results
    Tracker::VertexLilst vertex_found = vertex_finder.GetResults();
    std::vector<Tracker::Vertex *> vertices;
    for (const auto &vertex : vertex_found)
    {
        vertices.push_back(vertex.get());
    }

    //--------------------------------------------------------------------------------------
    // Run on multiple events

    int N = args["nevents"].as<int>();
    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    float duration = 0;

    // Generate 1000 Events first
    std::vector<std::vector<Tracker::Track *>> tracklist_all;
    for (auto i = 0; i < N; i++)
    {
        hits = genHitsTracks(N_TRACKS_PER_EVENT);
        auto track_finder_temp = Tracker::TrackFinder(hits, false);
        track_finder_temp.FindAll();

        track_found = track_finder_temp.GetResults();
        tracks.clear();
        for (auto &track : track_found)
        {
            tracks.push_back(track.release());
        }
        tracklist_all.push_back(tracks);
    }

    // Run LS fit, cpa chi2 function
    print(util::py::f("\n\nRun track/vertex reconstruction on {} event. Using Least square fit, with CPA chi2", N));
    start = std::chrono::high_resolution_clock::now();
    TH1D *h1 = new TH1D("h1", "Number of vertex", 10, 0, 10);
    TH1D *h2 = new TH1D("h2", "Residual_r", 100, 0, 80);
    TH1D *h3 = new TH1D("h3", "Vertex chi2", 100, 0, (N_TRACKS_PER_EVENT * 3 - 4) * 5);
    for (int i = 0; i < N; i++)
    {
        std::cout << "Event " << i << "\r";

        vertex_finder = Tracker::VertexFinder(tracklist_all[i], false);
        vertex_finder.FindAll();
        vertex_found = vertex_finder.GetResults();

        h1->Fill(vertex_found.size());

        for (const auto &vertex : vertex_found)
        {
            h2->Fill(vertex->params.segment(0, 3).norm());
            h3->Fill(vertex->chi2);
        }
    }
    std::cout << std::endl;
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(stop - start).count();
    std::cout << duration << " seconds for " << N << " events" << std::endl;

    // Plot and save output
    util::io::create_directory("plots");

    TCanvas *c1 = new TCanvas("c1", "N vertex", 800, 600);
    h1->SetLineColor(gStyle->GetColorPalette(0));
    h1->Draw();
    c1->BuildLegend();
    c1->Update();
    c1->SaveAs("plots/test4_nvertex.pdf");

    TCanvas *c2 = new TCanvas("c2", "Residual R", 800, 600);
    h2->SetLineColor(gStyle->GetColorPalette(0));
    h2->Draw();
    c2->BuildLegend();
    c2->Update();
    c2->SaveAs("plots/test4_residual_r.pdf");

    TCanvas *c3 = new TCanvas("c3", "chi2", 800, 600);
    h3->SetLineColor(gStyle->GetColorPalette(0));
    h3->Draw();
    c3->BuildLegend();
    c3->Update();
    c3->SaveAs("plots/test4_chi2.pdf");

    return 0;
}
