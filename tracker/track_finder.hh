#ifndef track_finder_HH
#define track_finder_HH

#include <iostream>
#include <memory>

#include "types.hh"
#include "util.hh"

namespace Tracker
{

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
        int FindOnce(TrackSeed *seed);

        // Find all tracks
        // Return number of tracks found
        int FindAll();
        TrackList GetResults() {return std::move(tracks_found);}

        // Remove used hits in hits list and seed list
        int RemoveUsed();

        // Get the event summary
        std::string Summary();

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
} // namespace Tracker

#endif