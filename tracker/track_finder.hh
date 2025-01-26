#ifndef track_finder_HH
#define track_finder_HH

#include <iostream>
#include <memory>

#include "types.hh"
#include "util.hh"

namespace Tracker
{
    class TrackSeed
    {
    public:
        TrackSeed(DigiHit* hit1, DigiHit* hit2) 
        {
            float c = 299.7; // speed of light [mm/ns]
            VectorXd dvec = hit2->vec4.array() - hit1->vec4.array();
            this->dstep = std::abs(hit2->get_step() - hit1->get_step());
            this->dr = dvec.segment(0,3).norm();
            this->dt = std::abs(dvec(3));
            this->score = std::abs(dr/c - dt);

            // Make a pair and put the earlier one in front
            this->hits = hit1->t() < hit2->t() ? std::make_pair(hit1, hit2) : std::make_pair(hit2, hit1);
        }

        // Required data
        std::pair<DigiHit*, DigiHit*> hits;
        float score;
        float dr,dt,dstep;

        // float GetScore() { return score; }
    };

    class TrackFinder
    {
    public:
        TrackFinder(std::vector<DigiHit *> allHits,
                    bool debug = false,
                    bool debug_kalman = false);

        // Clear the internal states
        void Clear();

        // Make track seeds
        void MakeSeeds();

        // Group hits by layer
        void GroupHits();

        // Find track for one seed 
        // Resutls are saved in a class variable `hits_found_temp`
        // Return the status:
        // 0: succeed
        // -1: not enough hits.
        // -2: track rejected by chi2 cut
        // -3: track rejected by speed cut
        int FindOnce(TrackSeed* seed);

        // Find all tracks
        // Return number of tracks found
        int FindAll();

        // Remove used hits in hits list and seed list
        int RemoveUsed();

    protected:
        // python-like print function, enable for debug mode only.
        template <typename... Args>
        inline void print_dbg(Args... args) {if (DEBUG) print(args...);}

        bool DEBUG, DEBUG_KALMAN;

        // All available hits and seeds
        std::vector<DigiHit *> hits_all;
        std::vector<TrackSeed *> seeds_unused;
        std::unordered_map<int, std::vector<DigiHit*>> hits_grouped;
        std::vector<int> hits_groups;

        // Temporary holders
        std::vector<DigiHit*> hits_found_temp;
        TrackList tracks_found;

        // Tracker configuration parameters
        std::unordered_map<std::string, double> config;
        int status; 
    };
} // namespace Tracker

#endif