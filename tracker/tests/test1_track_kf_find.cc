// stdlib
#include <cmath>
#include <chrono>

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>


#include "libs/cxxopts.hpp"
#include "test_util.hh"


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
    std::vector<Tracker::DigiHit *> hits = genHitsTracks(40);
    auto track_finder = Tracker::TrackFinder(false);
    track_finder.FindAll(hits);
    auto summary = track_finder.Summary();
    print(summary);

    int N = 100;
    VectorXd chi2(N);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++)
    {
        // std::cout << "--------------------------- " << i << "------------------" << std::endl;
        hits = genHitsTracks(40);
        auto track_finder = Tracker::TrackFinder(false);
        track_finder.FindAll(hits);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;

    // // std::cout << "chi2: " << chi2.transpose() << std::endl;
    // std::cout << "chi2 mean: " << chi2.mean() << std::endl;

    return 0;
}
