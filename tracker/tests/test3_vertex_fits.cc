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
#include "util_globals.hh"

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

    int N_TRACKS_PER_EVENT = args["ntracks"].as<int>();
    // clang-format on

    // Look at one track
    print("Run track/vertex reconstruction on one event.");
    std::vector<Tracker::DigiHit *> hits = genHitsTracks(N_TRACKS_PER_EVENT, {1000,1000,-8000, 0});
    auto track_finder = Tracker::TrackFinder(false);
    track_finder.FindAll(hits);
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

    // Run LS fit
    auto vertex_fitter = Kalman::LSVertex4DFitter(3);
    vertex_fitter.fit(tracks);
    print("Vertex LS Fit result (x,y,z,t):", vertex_fitter.params.transpose());
    print("Vertex LS Fit cov (x,y,z,t)\n", vertex_fitter.cov);

    // Run KF fit
    auto vertex_fitter_kf = std::make_unique<Kalman::KalmanVertex4D>(false);
    auto vertex = vertex_fitter_kf->run_filter(tracks);
    print("Vertex KF Fit result (x,y,z,t):", vertex->params.transpose());
    print("Vertex LS Fit cov (x,y,z,t)\n", vertex->cov);

    int N = args["nevents"].as<int>();
    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    float duration = 0;

    // // Generate 1000 Events first
    // std::vector<std::vector<Tracker::Track *>> tracklist_all;
    // for (auto i = 0; i < N; i++)
    // {
    //     hits = genHitsTracks(N_TRACKS_PER_EVENT, {1000,1000,-8000, 0});
    //     auto track_finder_temp = Tracker::TrackFinder(false);
    //     track_finder_temp.FindAll(hits);

    //     track_found = track_finder_temp.GetResults();
    //     tracks.clear();
    //     for (auto &track : track_found)
    //     {
    //         tracks.push_back(track.release());
    //     }
    //     tracklist_all.push_back(tracks);
    // }

    // // Run LS fit, cpa chi2 function
    // print(util::py::f("\n\nRun track/vertex reconstruction on {} event. Using Least square fit, with CPA chi2", N));
    // start = std::chrono::high_resolution_clock::now();
    // TH1D *h1 = new TH1D("h1", "Residual_r LS cpa", 100, 0, 80);
    // TH1D *h2 = new TH1D("h2", "Residual_t LS cpa", 100, -7, 7);
    // TH1D *h3 = new TH1D("h3", "chi2 LS cpa", 100, 0, (N_TRACKS_PER_EVENT * 3 - 4) * 5);
    // for (int i = 0; i < N; i++)
    // {
    //     std::cout << "Event " << i << "\r";

    //     vertex_fitter.fit(tracklist_all[i]);
    //     h1->Fill(vertex_fitter.params.segment(0, 3).norm());
    //     h2->Fill(vertex_fitter.params(3));
    //     h3->Fill(vertex_fitter.chi2);
    // }
    // std::cout << std::endl;
    // stop = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration<float>(stop - start).count();
    // std::cout << duration << " seconds for " << N << " events" << std::endl;

    // // Run LS fit, same-time chi2 function
    // auto vertex_fitter_time = Kalman::LSVertex4DFitter(2);
    // print(util::py::f("\n\nRun track/vertex reconstruction on {} event. Using Least square fit, with same-time chi2", N));
    // start = std::chrono::high_resolution_clock::now();
    // TH1D *h4 = new TH1D("h4", "Residual_r LS time", 100, 0, 80);
    // TH1D *h5 = new TH1D("h5", "Residual_t LS time", 100, -7, 7);
    // TH1D *h6 = new TH1D("h6", "chi2 LS time", 100, 0, (N_TRACKS_PER_EVENT * 3 - 4) * 5);
    // for (int i = 0; i < N; i++)
    // {
    //     std::cout << "Event " << i << "\r";

    //     vertex_fitter_time.fit(tracklist_all[i]);
    //     h4->Fill(vertex_fitter_time.params.segment(0, 3).norm());
    //     h5->Fill(vertex_fitter_time.params(3));
    //     h6->Fill(vertex_fitter_time.chi2);
    // }
    // std::cout << std::endl;
    // stop = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration<float>(stop - start).count();
    // std::cout << duration << " seconds for " << N << " events" << std::endl;

    // // Run KF fit
    // print(util::py::f("\n\nRun track/vertex reconstruction on {} event. Using Kalman Filter fit.", N));
    // start = std::chrono::high_resolution_clock::now();
    // TH1D *h7 = new TH1D("h7", "Residual_r KF", 100, 0, 80);
    // TH1D *h8 = new TH1D("h8", "Residual_t KF", 100, -7, 7);
    // TH1D *h9 = new TH1D("h9", "chi2 KF", 100, 0, (N_TRACKS_PER_EVENT * 3 - 4) * 5);
    // for (int i = 0; i < N; i++)
    // {
    //     std::cout << "Event " << i << "\r";

    //     auto _vertex = vertex_fitter_kf->run_filter(tracklist_all[i]);
    //     h7->Fill(_vertex->params.segment(0, 3).norm());
    //     h8->Fill(_vertex->params(3));
    //     h9->Fill(_vertex->chi2);
    // }
    // std::cout << std::endl;
    // stop = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration<float>(stop - start).count();
    // std::cout << duration << " seconds for " << N << " events" << std::endl;

    // print("\n");

    // gROOT->SetStyle("Modern");
    // // gStyle->SetPalette(kRainBow);
    // setTableauPalette();


    // // Plot and save output
    // util::io::create_directory("plots");

    // TCanvas *c1 = new TCanvas("c1", "Residual R", 800, 600);
    // h1->SetLineColor(gStyle->GetColorPalette(0));
    // h1->Draw();
    // h4->SetLineColor(gStyle->GetColorPalette(1));
    // h4->Draw("SAME");
    // h7->SetLineColor(gStyle->GetColorPalette(2));
    // h7->Draw("SAME");
    // c1->BuildLegend();
    // c1->Update();
    // c1->SaveAs("plots/test3_residual_r.pdf");

    // TCanvas *c2 = new TCanvas("c2", "Residual t", 800, 600);
    // h2->SetLineColor(gStyle->GetColorPalette(0));
    // h2->Draw();
    // h5->SetLineColor(gStyle->GetColorPalette(1));
    // h5->Draw("SAME");
    // h8->SetLineColor(gStyle->GetColorPalette(2));
    // h8->Draw("SAME");
    // c2->BuildLegend();
    // c2->Update();
    // c2->SaveAs("plots/test3_residual_t.pdf");

    // TCanvas *c3 = new TCanvas("c3", "chi2", 800, 600);
    // h3->SetLineColor(gStyle->GetColorPalette(0));
    // h3->Draw();
    // h6->SetLineColor(gStyle->GetColorPalette(1));
    // h6->Draw("SAME");
    // h9->SetLineColor(gStyle->GetColorPalette(2));
    // h9->Draw("SAME");
    // c3->BuildLegend();
    // c3->Update();
    // c3->SaveAs("plots/test3_chi2.pdf");

    return 0;
}
