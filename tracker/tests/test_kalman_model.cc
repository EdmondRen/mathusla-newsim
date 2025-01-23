#include <cmath>

#include "../kalman_model.hh"

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

TRandom3 rng;

Tracker::HitList genHits(double x0_span, double y0_span, double z0, double Ax_span, double Ay_span, double v = 30, int n_layers = 6)
{
    // Generate a vector with 10 elements, linearly spaced from 0 to 9
    auto layer_ids = Eigen::VectorXd::LinSpaced(n_layers, 0, n_layers - 1);
    auto layer_zs = layer_ids * 80; // y_coordinates

    float det_width = 3.5;
    float det_height = 1;
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

    Tracker::HitList hits;
    for (int i = 0; i < n_layers; i++)
    {
        auto dz = layer_zs[i] - layer_zs[0];
        double hit_truth[6] = {x0 + Ax * dz, y0 + Ay * dz, t0 + At * dz, Ax, Ay, At};
        if (i % 2 == 0)
            hits.push_back(new Tracker::DigiHit(round(hit_truth[0] / det_width) * det_width,
                                                hit_truth[1] + rng.Gaus(0, unc_longi),
                                                layer_zs[i],
                                                hit_truth[2], unc_trans, unc_longi, unc_verti, unc_time,
                                                independent_var_ind, 0));
        else
            hits.push_back(new Tracker::DigiHit(hit_truth[0] + rng.Gaus(0, unc_longi),
                                                round(hit_truth[1] / det_width) * det_width,
                                                layer_zs[i],
                                                hit_truth[2], unc_longi, unc_trans, unc_verti, unc_time,
                                                independent_var_ind, 0));
    }
    return hits;
}

int main()
{

    rng.SetSeed(0);

    // Look at one track
    Tracker::HitList hits = genHits(30, 30, 0, 0.1, 0.1);
    auto fitter = new Kalman::KalmanTrack4D(0, 4, 6, true);
    auto track = fitter->run_filter(hits);

    hits = genHits(30, 30, 0, 0.1, 0.1);
    fitter->run_filter(hits);

    int N = 10;
    VectorXd chi2(N);
    for (int i = 0; i < N; i++)
    {   
    std::cout << "--------------------------- " << i << "------------------" << std::endl;
        // auto fitter_new = Kalman::KalmanTrack4D(0, 4, 6, true);
        auto hits_temp = genHits(30, 30, 0, 0.1, 0.1);
        auto track2 = fitter->run_filter(hits_temp);
        chi2(i) = track2->chi2;
        delete track2;
        for (auto hit:hits_temp)
            delete hit;
    }

    std::cout << "chi2: " << chi2.transpose() << std::endl;
    std::cout << "chi2 mean: " << chi2.mean() << std::endl;

    return 0;
}