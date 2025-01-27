// stdlib
#include <cmath>
#include <chrono>

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

#include "util.hh"
#include "libs/cxxopts.hpp"

#include "../track_finder.hh"
#include "../kalman_model.hh"

TRandom3 rng;

std::vector<Tracker::DigiHit *> genHits(double x0_span, double y0_span, double z0, double Ax_span, double Ay_span, double v = 300, int n_layers = 6)
{
    auto layer_ids = Eigen::VectorXd::LinSpaced(n_layers, 0, n_layers - 1);
    auto layer_gap = 800.;                                // 80 cm
    auto layer_zs = (layer_ids * layer_gap).array() + z0; // y_coordinates

    float det_width = 35;
    float det_height = 10;
    float refraction_index = 1.89;

    float unc_time = 1; // ns, coincidence resolution
    float unc_trans = det_width / std::sqrt(12.);
    float unc_longi = unc_time * v / refraction_index;
    float unc_verti = det_height / std::sqrt(12.);

    auto x0 = rng.Uniform(x0_span * 2) - x0_span;
    auto y0 = rng.Uniform(y0_span * 2) - y0_span;
    auto t0 = 0.0;
    auto Ax = rng.Uniform(Ax_span * 2) - Ax_span;
    auto Ay = rng.Uniform(Ay_span * 2) - Ay_span;
    auto At = std::sqrt(1. + std::pow(Ax, 2) + std::pow(Ay, 2)) / v;

    int independent_var_ind = 2;

    std::vector<Tracker::DigiHit *> hits;
    for (int i = 0; i < n_layers; i++)
    {
        auto dz = layer_zs[i];
        double hit_truth[6] = {x0 + Ax * dz, y0 + Ay * dz, t0 + At * dz, Ax, Ay, At};
        Tracker::DigiHit *hit;
        if (i % 2 == 0)
            hit = new Tracker::DigiHit(round(hit_truth[0] / det_width) * det_width,
                                       hit_truth[1] + rng.Gaus(0, unc_longi),
                                       layer_zs[i],
                                       hit_truth[2] + rng.Gaus(0, unc_time),
                                       unc_trans, unc_longi, unc_verti, unc_time,
                                       independent_var_ind, i, 0);
        else
            hit = new Tracker::DigiHit(hit_truth[0] + rng.Gaus(0, unc_longi),
                                       round(hit_truth[1] / det_width) * det_width,
                                       layer_zs[i],
                                       hit_truth[2] + rng.Gaus(0, unc_time),
                                       unc_longi, unc_trans, unc_verti, unc_time,
                                       independent_var_ind, i, 0);
        hit->id = i;
        hits.push_back(hit);
    }
    return hits;
}

std::vector<Tracker::DigiHit *> genHitsTracks(int ntracks)
{
    std::vector<Tracker::DigiHit *> hits;

    for (int i = 0; i < ntracks; i++)
    {
        auto hits_temp = genHits(0, 0, 4000, 1, 1);
        hits.insert(hits.end(), std::make_move_iterator(hits_temp.begin()), std::make_move_iterator(hits_temp.end()));
    }

    for (size_t i = 0; i < hits.size(); i++)
    {
        hits[i]->id = i;
    }
    return hits;
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
    // clang-format on    

    rng.SetSeed(args["seed"].as<int>());

    // Look at one track
    // std::vector<Tracker::DigiHit *> hits = genHits(300, 300, 0, 0.3, 0.3);
    std::vector<Tracker::DigiHit *> hits = genHitsTracks(10);
    auto track_finder = Tracker::TrackFinder(hits, false);
    track_finder.FindAll();
    auto summary = track_finder.Summary();
    print(summary);

    // Make a list of track pointers
    Tracker::TrackList track_found = track_finder.GetResults();
    std::vector<Tracker::Track *> tracks;
    for (const auto & track: track_found)
    {
        track->convert_to_time();
        tracks.push_back(track.get());
        }



    auto vertex_fitter = Kalman::LSVertex4DFitter();  
    vertex_fitter.fit(tracks);

    print("Vertex Fit result (x,y,z,t)", vertex_fitter.parameters); 

    int N = 100;
    VectorXd chi2(N);

    // auto start = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < N; i++)
    // {
    //     // std::cout << "--------------------------- " << i << "------------------" << std::endl;
    //     hits = genHitsTracks(40);
    //     auto track_finder = Tracker::TrackFinder(hits, false);
    //     track_finder.FindAll();

    //      track_found.clear();
    //     tracks.clear();
    //     for (const auto & track: track_found)
    //     {
    //         track->convert_to_time();
    //         tracks.push_back(track.get());
    //     }
    //     vertex_fitter.fit(tracks);

    // }
    // auto stop = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << duration.count() << std::endl;

    // // std::cout << "chi2: " << chi2.transpose() << std::endl;
    // std::cout << "chi2 mean: " << chi2.mean() << std::endl;

    return 0;
}
