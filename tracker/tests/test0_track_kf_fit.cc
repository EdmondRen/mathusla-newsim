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


void run_once()
{
    std::vector<Tracker::DigiHit*> hits = genHits(30, 30, 0, 0.1, 0.1);
    auto fitter = std::make_unique<Kalman::KalmanTrack4D>(0, 4, 6, false);
    auto track = fitter->run_filter(hits);
    auto chi2 = track->chi2;
    std::cout << "chi2: " << chi2 << std::endl;
}

int main()
{

    rng.SetSeed(0);

    // Look at one track
    std::vector<Tracker::DigiHit*> hits = genHits(30, 30, 0, 0.1, 0.1);
    auto fitter = std::make_unique<Kalman::KalmanTrack4D>(0, 4, 6, false);
    auto track = fitter->run_filter(hits);

    // hits = genHits(30, 30, 0, 0.1, 0.1);
    // fitter->run_filter(hits);



    int N = 1000;
    auto start = std::chrono::high_resolution_clock::now();
    VectorXd chi2(N);

    for (int i = 0; i < N; i++)
    {
        // std::cout << "--------------------------- " << i << "------------------" << std::endl;
        // auto fitter_new = Kalman::KalmanTrack4D(0, 4, 6, true);
        auto hits_temp = genHits(30, 30, 0, 0.1, 0.1);
        auto track2 = fitter->run_filter(hits_temp);
        chi2(i) = track2->chi2;
        // delete track2;
        // for (auto hit : hits)
        //     std::cout << "hitx: " << hit->x() << std::endl;
        // run_once();
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;


    // std::cout << "chi2: " << chi2.transpose() << std::endl;
    std::cout << "chi2 mean: " << chi2.mean() << std::endl;

    return 0;
}
