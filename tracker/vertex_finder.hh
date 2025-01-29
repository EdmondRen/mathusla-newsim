#ifndef vertex_finder_HH
#define vertex_finder_HH

#include <iostream>
#include <memory>

#include "types_and_helper.hh"
#include "util.hh"
#include "kalman_model.hh"

namespace Tracker
{
    using Kalman::VertexSeed;

    class VertexFinder
    {
    public:
        VertexFinder(std::vector<Track *> allTracks,
                    bool debug = false,
                    bool debug_kalman = false);

        // Clear the internal states
        void Clear();

        // Make vertex seeds
        void MakeSeeds();

        // Find vertex for one seed
        // Resutls are saved in a class variable `hits_found_temp`
        // Return the status:
        // 0: succeed
        // -1: not enough hits.
        // -2: track rejected by chi2 cut
        // -3: track rejected by speed cut
        int FindOnce(VertexSeed *seed);

        // Find all vertices
        // Return number of vertices found
        int FindAll();
        VertexLilst GetResults() {return std::move(vertices_found);}

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
        std::vector<Track *> tracks_all;
        std::unordered_map<int, Track *> tracks_unused;
        std::vector<Kalman::VertexSeed *> seeds_unused;

        // Temporary holders
        std::vector<Track *> tracks_found_temp;
        VertexLilst vertices_found;

        // configuration parameters
        std::unordered_map<std::string, double> config;
        int status;

        // Summary of vertex finding results
        std::string summary;
    }; // class VertexFinder
} // namespace Tracker

#endif