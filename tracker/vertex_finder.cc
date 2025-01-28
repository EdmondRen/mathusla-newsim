// ROOT libs
#include <TMath.h>
#include <Math/ProbFunc.h>

#include "vertex_finder.hh"
#include "kalman_model.hh"

namespace Tracker
{

    VertexFinder::VertexFinder(std::vector<Track *> allTracks,
                               bool debug, bool debug_kalman) : DEBUG(debug), DEBUG_KALMAN(debug_kalman), tracks_all(allTracks)
    {
        this->tracks_found_temp = std::vector<Track *>();
        this->vertices_found = VertexLilst();

        // Set the default parameters
        config["vertex_cut_SeedDist"] = 300;                               // [mm] Maxmum distance between two tracks
        config["vertex_cut_SeedChi2"] = 50;                                // maximum chi2 of the seed
        config["vertex_cut_TrackAddDist"] = 300;                           //  Distance cut to add a track
        config["vertex_cut_TrackAddChi2"] = 9;                             // chi2 cut to add a track
        config["vertex_cut_TrackDropChi2"] = 3;                            // chi2 cut to drop a track
        config["vertex_cut_VertexChi2Reduced"] = 5;                        // Vertex chi2-square cut
        config["vertex_cut_VertexNHitsMin"] = 2;                           // Vertex minimum number of tracks
        config["vertex_fit_MultipleScattering"] = 0;                       // Account for multiple scattering in fit. 0: OFF
        config["vertex_multiple_scattering_p"] = 500;                      // Particle momentum asumed for multiple scattering
        config["vertex_multiple_scattering_length"] = 0.06823501107481977; // Material thickness in the unit of scattering length. Calculated for 1 cm plastic scintillator and 0.6 mm aluminum.
    }

    void VertexFinder::Clear()
    {
    }

    void VertexFinder::MakeSeeds()
    {
        for (size_t i = 0; i < tracks_all.size() - 1; i++)
        {
            for (size_t j = i + 1; j < tracks_all.size(); j++)
            {

                auto new_seed = new Kalman::VertexSeed(tracks_all[i], tracks_all[j]);
                if (new_seed->distance < config["vertex_cut_SeedDist"])
                {
                    if (new_seed->GetChi2() < config["vertex_cut_SeedChi2"])
                        this->seeds_unused.push_back(new_seed);
                }
                else
                    delete new_seed;
            }
        }

        // Calculate score of the seed
        for (auto seed : seeds_unused)
        {
            // Check how many tracks could be associated with each seed
            for (auto track : tracks_all)
            {
                if (track->get_same_time_dist(seed->vertex_fit) < config["vertex_cut_TrackAddDist"])
                    seed->ntracks_found += 1;
            }
        }

        // Sort seeds by score from small to large. Larger score is better
        // Put beter ones at the end
        std::sort(seeds_unused.begin(), seeds_unused.end(),
                  [](Kalman::VertexSeed *a, Kalman::VertexSeed *b) -> bool
                  { return a->score < b->score; });

        int iprint = 0;
        print_dbg(util::py::f("Making seeds, total {} seeds. Printing best 30. Smaller score is better.", seeds_unused.size()));
        if (DEBUG)
        {
            for (int i = seeds_unused.size() - 1; i >= 0; --i)
            {
                iprint += 1;
                auto seed = seeds_unused[i];
                print_dbg(util::py::f("  Seed [{:03d}, {:03d}]: score {:.4f}, ntracks_found: {:.0f}, distance {:.3f}, chi2 {:.3f}",
                                      seed->tracks.first->id, seed->tracks.second->id, seed->score, seed->ntracks_found, seed->distance, seed->chi2));
                if (iprint > 30)
                    break;
            }
        }

        // Make a map of all tracks
        for (auto track : tracks_all)
            tracks_unused[track->id] = track;
    }

    int VertexFinder::FindOnce(VertexSeed *seed)
    {
        print_dbg(util::py::f("\nUsing Seed [{},{}]", seed->tracks.first->id, seed->tracks.second->id));

        tracks_found_temp.clear();
        auto finder = Kalman::KalmanVertex4D(false); // config["track_fit_MultipleScattering"]

        // Initialize the finder with seeds
        tracks_found_temp.push_back(seed->tracks.first);
        tracks_found_temp.push_back(seed->tracks.second);
        finder.init_state(*seed);

        // Calculate the chi2 contribution of each remaining track
        // Only keep the minimum one
        int chi2_min_ind = -1;
        float chi2_min_val = std::numeric_limits<float>::infinity();

        for (const auto &pair : tracks_unused)
        {
            auto track_id = pair.first;
            auto track = pair.second;
            

            float track_chi2;
            track_chi2 = finder.try_measurement(*track,
                                                config["vertex_cut_TrackAddDist"],
                                                config["vertex_cut_TrackAddChi2"]);
            // print_dbg("    hit chi2:", hit_chi2);
            // find the smallest one
            if (track_chi2 >= 0 && track_chi2 < chi2_min_val)
            {
                chi2_min_ind = track_id;
                chi2_min_val = track_chi2;
            }
        }
        print_dbg(util::py::f("  Minimum chi2 is {:.3f}, from track {}.", chi2_min_val, chi2_min_ind));

        // Add the hit with minimum chi2 and update the finder status
        if (chi2_min_ind != -1)
        {
            tracks_found_temp.push_back(tracks_unused[chi2_min_ind]); // add to list
            finder.add_measurement(*tracks_unused[chi2_min_ind]);     // update finder
        }

        return tracks_found_temp.size();
    }

    int VertexFinder::FindAll()
    {
        // Sort hits into groups
        print_dbg(util::py::f("Event contains {} tracks", tracks_all.size()));

        // Seeding
        MakeSeeds();

        if (tracks_unused.size() < config["vertex_cut_VertexNHitsMin"])
            return 0;

        // Loop though required number of hits
        this->tracks_found_temp.clear();
        int ntracks_min = config["vertex_cut_VertexNHitsMin"];
        while (this->seeds_unused.size() > 0)
        {
            // Round 1: find hits based on seed
            int ntracks_found = FindOnce(seeds_unused.back());

            // Cut on number of tracks
            if (ntracks_found < ntracks_min)
            {
                print_dbg(util::py::f("-x Vertex rejected, only {} tracks", ntracks_found));
                seeds_unused.pop_back();
                continue;
            }

            std::string tracks_found_ids = "";
            for (auto track : tracks_found_temp)
            {
                tracks_found_ids += std::to_string(track->id) + " ";
            }
            print_dbg(util::py::f("-> Found vertex with {} tracks:", tracks_found_ids.size()), tracks_found_ids);

            // Round 2: Drop outliers
            // Discard this track if the seed contributes too much to chi2

            // Round 3: Run fit again using least square minimizer
            auto vertex_fitter = Kalman::LSVertex4DFitter(1);
            auto vertex_found = vertex_fitter.run_fit(tracks_found_temp);
            print_dbg("Vertex LS Fit result (x,y,z,t):", vertex_found->params.transpose());

            // Cut on chi2
            int ndof = tracks_found_temp.size() * 3 - 4;
            float chi2 = vertex_found->chi2;
            float chi2_prob = ROOT::Math::chisquared_cdf(chi2, ndof);
            bool passed = ((chi2 / ndof) < config["vertex_cut_VertexChi2Reduced"]);
            if (!passed)
            {
                print_dbg(util::py::f("-x Vertex rejected, chi2/ndof ({}/{}) exceeds limit, prob = {}", chi2, ndof, 1 - chi2_prob));
                seeds_unused.pop_back();
                continue;
            }

            // Setting the metadata of this vertex
            vertex_found->id = this->vertices_found.size();
            vertex_found->track_ids = std::vector<int>();

            tracks_found_ids = "";
            for (auto hit : tracks_found_temp)
            {
                tracks_found_ids += std::to_string(hit->id) + " ";
                vertex_found->track_ids.push_back(hit->id);
            }
            print_dbg(util::py::f("-> Found vertex with {} tracks:", tracks_found_temp.size()), tracks_found_ids);

            // Add vertex to the list
            this->vertices_found.push_back(std::move(vertex_found));
            print_dbg(util::py::f("-> Vertex #{} added: chi2 {:.3f}", this->vertices_found.size(), this->vertices_found.back()->chi2));

            // Finally, remove other seeds that have hits in this vertex
            RemoveUsed();
            print_dbg("  Remaining seeds:", seeds_unused.size());
        }

        return -1;
    }

    int VertexFinder::RemoveUsed()
    {
        // removing used seeds
        std::vector<VertexSeed *> seeds_unused_temp;
        std::string seeds_removed_list;
        for (size_t i = 0; i < seeds_unused.size(); i++)
        {
            bool veto_this = false;
            auto seed = seeds_unused[i];
            for (auto track : tracks_found_temp)
            {
                if (seed->tracks.first->id == track->id ||
                    seed->tracks.second->id == track->id)
                {
                    veto_this = true;
                    break;
                }
            }
            if (veto_this)
            {
                seeds_removed_list += "[" + std::to_string(seed->tracks.first->id) + "," + std::to_string(seed->tracks.second->id) + "] ";
                delete seed;
            }
            else
                seeds_unused_temp.push_back(seed);
        }
        print_dbg(util::py::f("  Seed removed: {}", seeds_removed_list));
        this->seeds_unused = seeds_unused_temp;

        // removing used tracks
        for (auto track : tracks_found_temp)
        {
            tracks_unused.erase(track->id);
        }

        return 0;
    }

    std::string VertexFinder::Summary()
    {

        summary = "";
        summary += util::py::f("Event vertex-finding summary:\n");
        summary += util::py::f("  Total tracks: {}\n", tracks_all.size());
        summary += util::py::f("  Total vertices: {}\n", vertices_found.size());
        for (auto &vertex : vertices_found)
        {
            summary += util::py::f("   vertex: with {} hits\n", vertex->track_ids.size());
        }

        return summary;
    }

} // namespace Tracker
