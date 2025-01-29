#ifndef test_util_hh
#define test_util_hh

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

#include "../track_finder.hh"
#include "../kalman_model.hh"

#include "util.hh"


extern TRandom3 rng;

std::vector<Tracker::DigiHit *> genHits(double x0_span, double y0_span, double z0, double Ax_span, double Ay_span, double v = 300, int n_layers = 6, const std::vector<double> &vertex_pos = {0,0,0,0});

std::vector<Tracker::DigiHit *> genHitsTracks(int ntracks, const std::vector<double> &vertex_pos =  std::vector<double>{0,0,0,0});

// std::vector<Tracker::DigiHit *> genHitsVertices(int ntracks);



#endif
