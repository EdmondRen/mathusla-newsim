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

    auto x0 = rng.Uniform(x0_span);
    auto y0 = rng.Uniform(y0_span);
    auto t0 = 0.0;
    auto Ax = rng.Uniform(Ax_span);
    auto Ay = rng.Uniform(Ay_span);
    auto At = std::sqrt(1. + std::pow(Ax, 2) + std::pow(Ay, 2)) / 2;

    int independent_var_ind = 2;

    Tracker::HitList hits;
    for (int i = 0; i < n_layers; i++)
    {
        auto dz = layer_zs[i] - layer_zs[0];
        float hit_truth[6] = {x0, y0, t0, Ax, Ay, At};
        if (i % 2 == 0)
            hits.push_back(new Tracker::DigiHit(round(x0 / det_width) * det_width,
                                                y0 + rng.Gaus(0, unc_longi),
                                                layer_zs[i],
                                                t0, unc_trans, unc_longi, unc_verti, unc_time,
                                                independent_var_ind, 0));
        else
            hits.push_back(new Tracker::DigiHit(x0 + rng.Gaus(0, unc_longi),
                                                round(y0 / det_width) * det_width,
                                                layer_zs[i],
                                                t0, unc_trans, unc_longi, unc_verti, unc_time,
                                                independent_var_ind, 0));
    }
    return hits;
}

int main()
{

    rng.SetSeed(0);

    auto hits = genHits(30, 30, 0, 0.1, 0.1);

    auto fitter = Kalman::KalmanTrack4D(0);
    auto track = fitter.run_filter(hits);

    return 0;
}