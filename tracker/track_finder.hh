#ifndef track_finder_HH
#define track_finder_HH

#include <iostream>
#include <memory>

#include "kalman_model.hh"
#include "types_and_helper.hh"
#include "util.hh"

namespace Tracker
{

    class TrackFinder
    {
    public:
        TrackFinder(bool debug = false,
                    bool debug_kalman = false);

        // Clear the internal states
        void Clear();

        // Make track seeds
        void MakeSeeds();
        void SortSeedsOccurance();

        // Group hits by layer
        void GroupHits();

        // Find track for one seed
        // Resutls are saved in a class variable `hits_found_temp`
        // Return the status:
        // 0: succeed
        // -1: not enough hits.
        // -2: track rejected by chi2 cut
        // -3: track rejected by speed cut
        int FindOnce(TrackSeed *seed);

        // Find all tracks
        // Return number of tracks found
        int FindAll(std::vector<DigiHit *> allHits);
        TrackList GetResults() { return std::move(tracks_found); }

        // Remove used hits in hits list and seed list
        int RemoveUsed();

        // Get the event summary
        std::string Summary();
        int info_ntracks;
        int info_nhits;

        // Configuration
        void Config(std::map<std::string, double> &config_ext);

        double hits_found_temp_chi2;

    protected:
        // python-like print function, enable for debug mode only.
        template <typename... Args>
        inline void print_dbg(Args... args)
        {
            if (DEBUG)
                print(args...);
        }

        bool DEBUG, DEBUG_KALMAN;

        // All available hits and seeds
        std::vector<DigiHit *> hits_all;
        std::vector<TrackSeed *> seeds_unused;
        std::vector<TrackSeed *> seeds_unused_next;
        std::unordered_map<int, std::unordered_map<int, DigiHit *>> hits_grouped;
        std::vector<int> hits_groups;

        // Temporary holders
        std::vector<DigiHit *> hits_found_temp;
        TrackList tracks_found;

        // Tracker configuration parameters
        std::unordered_map<std::string, double> config;
        int status;

        // Summary of track finding results
        std::string summary;
    };


    class VertexTrackSeed
    {
    public:
        VertexTrackSeed(Vector4d point, DigiHit *hit2) : nhits_found(-1)
        {
            float c = 299.7; // speed of light [mm/ns]
            this->dvec = hit2->vec4 - point;
            this->dstep = std::abs(hit2->get_step() - point(hit2->iv_index));
            this->dr = dvec.segment(0, 3).norm();
            this->dt = std::abs(dvec(3));
            this->score = std::abs(dr / c - dt);

            // Make a pair and put the earlier one in front
            this->hits = std::make_pair(hit2->id, hit2);

            this->chi2prob_found = 1; // Default value of chi2 prob if not specified.
        }

        // Required data
        std::pair<int, DigiHit *> hits;
        float score;
        VectorXd dvec;
        float dr, dt, dstep;
        int nhits_found; // Number of hits found using this seed. Save time when reusing this seed.
        double chi2prob_found; //
        int nhit_occur;
        int layer;

        // float GetScore() { return score; }
    };

    class VertexTrackFinder
    {
    public:
        VertexTrackFinder(bool debug = false,
                    bool debug_kalman = false);

        // Clear the internal states
        void Clear();

        // Make track seeds
        void MakeSeeds();
        void SortSeedsOccurance();

        // Group hits by layer
        void GroupHits();

        // Find track for one seed
        // Resutls are saved in a class variable `hits_found_temp`
        // Return the status:
        // 0: succeed
        // -1: not enough hits.
        // -2: track rejected by chi2 cut
        // -3: track rejected by speed cut
        int FindOnce(VertexTrackSeed *seed);

        // Find all tracks
        // Return number of tracks found
        int FindAll(std::vector<DigiHit *> allHits, Vertex *_vertex_);
        TrackList GetResults() { return std::move(tracks_found); }

        // Remove used hits in hits list and seed list
        int RemoveUsed();

        // Get the event summary
        std::string Summary();
        int info_ntracks;
        int info_nhits;

        // Configuration
        void Config(std::map<std::string, double> &config_ext);

        double hits_found_temp_chi2;

    protected:
        bool DEBUG, DEBUG_KALMAN;
        Kalman::KalmanTrack4D *finder;

        // All available hits and seeds
        Vertex *vertex;
        std::vector<DigiHit *> hits_all;
        std::vector<VertexTrackSeed *> seeds_unused;
        std::vector<VertexTrackSeed *> seeds_unused_next;
        std::unordered_map<int, std::unordered_map<int, DigiHit *>> hits_grouped;
        std::vector<int> hits_groups;

        // Temporary holders
        std::vector<DigiHit *> hits_found_temp;
        TrackList tracks_found;

        // Tracker configuration parameters
        std::unordered_map<std::string, double> config;
        int status;

        // Summary of track finding results
        std::string summary;

        // python-like print function, enable for debug mode only.
        template <typename... Args>
        inline void print_dbg(Args... args)
        {
            if (DEBUG)
                print(args...);
        }        
    };    
} // namespace Tracker

#endif