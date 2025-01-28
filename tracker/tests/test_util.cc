#include "test_util.hh"


TRandom3 rng;

std::vector<Tracker::DigiHit *> genHits(double x0_span, double y0_span, double z0, double Ax_span, double Ay_span, double v, int n_layers)
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