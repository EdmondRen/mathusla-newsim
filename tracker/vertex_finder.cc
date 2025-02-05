// ROOT libs
#include <TMath.h>
#include <Math/ProbFunc.h>

#include "vertex_finder.hh"
#include "kalman_model.hh"

namespace Tracker
{

    VertexFinder::VertexFinder(bool debug, bool debug_kalman) : DEBUG(debug), DEBUG_KALMAN(debug_kalman)
    {
        this->tracks_found_temp = std::vector<Track *>();
        this->vertices_found = VertexLilst();

        // Set the default parameters
        config["vertex_cut_SeedDist"] = 3000;                              // [mm] Maxmum distance between two tracks
        config["vertex_cut_SeedChi2"] = 50;                                // maximum chi2 of the seed
        config["vertex_cut_TrackAddDist"] = 3000;                          //  Distance cut to add a track
        config["vertex_cut_TrackAddChi2"] = 15;                             // chi2 cut to add a track
        config["vertex_cut_TrackDropChi2"] = 7;                            // chi2 cut to drop a track
        config["vertex_cut_VertexChi2Reduced"] = 5;                        // Vertex chi2-square cut
        config["vertex_cut_VertexNHitsMin"] = 2;                           // Vertex minimum number of tracks
        config["vertex_fit_MultipleScattering"] = 1;                       // Account for multiple scattering in fit. 0: OFF
        config["vertex_multiple_scattering_p"] = 500;                      // Particle momentum asumed for multiple scattering
        config["vertex_multiple_scattering_length"] = 0.06823501107481977; // Material thickness in the unit of scattering length. Calculated for 1 cm plastic scintillator and 0.6 mm aluminum.
    }

    void VertexFinder::Config(std::map<std::string, double> &config_ext)
    {
        for (auto &pair : config)
        {
            if (config_ext.count(pair.first) > 0)
                config[pair.first] = config_ext[pair.first];
        }
    }

    void VertexFinder::Clear()
    {
    }

    void VertexFinder::MakeSeeds()
    {
        this->seeds_unused.clear();

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
        for (auto &seed : seeds_unused)
        {
            // Check how many tracks could be associated with each seed
            // Save the distance between each track and the seed
            for (auto &track : tracks_all)
            {
                if (track->id == seed->tracks.first->id ||
                    track->id == seed->tracks.second->id)
                    continue;

                double track_to_seed_dist = track->get_same_time_dist(seed->vertex_fit);
                seed->track_distance_list.push_back(std::make_pair(track->id, track_to_seed_dist));
                if (track_to_seed_dist < config["vertex_cut_TrackAddDist"])
                    seed->ntracks_found += 1;

                // Sort distance descending. Minimum distance at the end. Use backwards.
                std::sort(seed->track_distance_list.begin(), seed->track_distance_list.end(),
                          [](std::pair<int, double> &a, std::pair<int, double> &b)
                          { return a.second > b.second; });
            }
            

            // Score: defined as ntracks
            seed->score = seed->ntracks_found;
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

        // Try all tracks that are recorded in the seed
        for (int i = seed->track_distance_list.size() - 1; i >= 0; i--)
        {
            auto track_pair_dist = seed->track_distance_list[i];
            auto track_id = track_pair_dist.first;
            auto track_dist = track_pair_dist.second;

            if (track_dist > config["vertex_cut_TrackAddDist"] * 2 ||
                tracks_unused.count(track_id) == 0)
            {
                if (tracks_unused.count(track_id) != 0)
                {
                    print_dbg(util::py::f("  Track {} is ignored due to large distance of {} [mm]", track_id, track_dist));
                }

                // Delete it from the seed and skip this round of loop
                if (i != 0)
                {
                    seed->track_distance_list.erase(seed->track_distance_list.begin() + i);
                    continue;
                }
                else
                    break;
            }

            auto &track = tracks_unused[track_id];

            float track_chi2;
            track_chi2 = finder.try_measurement(*track,
                                                config["vertex_cut_TrackAddDist"],
                                                config["vertex_cut_TrackAddChi2"]);
            if (track_chi2 != -1)
            {
                print_dbg(util::py::f("  Track {} is added with chi2 {:.3f}.", track_id, track_chi2));
                tracks_found_temp.push_back(tracks_unused[track_id]); // add to list
                finder.add_measurement(*tracks_unused[track_id]);     // update finder
            }
        }

        return tracks_found_temp.size();
    }

    int VertexFinder::FindOnceLS(VertexSeed *seed)
    {
        print_dbg(util::py::f("\nUsing Seed [{},{}]", seed->tracks.first->id, seed->tracks.second->id));

        tracks_found_temp.clear();

        // Initialize the finder with seeds
        tracks_found_temp.push_back(seed->tracks.first);
        tracks_found_temp.push_back(seed->tracks.second);

        auto vertex_fitter = Kalman::LSVertex4DFitter(3, config["vertex_fit_MultipleScattering"]>0);
        auto vertex_found = vertex_fitter.run_fit(tracks_found_temp);
        float chi2_prev = vertex_found->chi2;

        // Try all tracks that are recorded in the seed
        for (int i = seed->track_distance_list.size() - 1; i >= 0; i--)
        {
            auto track_pair_dist = seed->track_distance_list[i];
            auto track_id = track_pair_dist.first;
            auto track_dist = track_pair_dist.second;

            if (track_dist > config["vertex_cut_TrackAddDist"] * 2 ||
                tracks_unused.count(track_id) == 0)
            {
                if (tracks_unused.count(track_id) != 0)
                {
                    print_dbg(util::py::f("  Track {} is ignored due to large distance of {} [mm]", track_id, track_dist));
                }

                // Delete it from the seed and skip this round of loop
                if (i != 0)
                {
                    seed->track_distance_list.erase(seed->track_distance_list.begin() + i);
                    continue;
                }
                else
                    break;
            }

            auto &track = tracks_unused[track_id];


            tracks_found_temp.push_back(tracks_unused[track_id]);
            auto vertex_temp = vertex_fitter.run_fit(tracks_found_temp);
            float track_chi2 = vertex_temp->chi2 - chi2_prev;

            if (track_chi2 < config["vertex_cut_TrackAddChi2"])
            {
                print_dbg(util::py::f("  Track {} is added with chi2 {:.3f}.", track_id, track_chi2));
                vertex_found = vertex_temp;
                chi2_prev = vertex_temp->chi2;
            }
            else
            {
                tracks_found_temp.pop_back(); // remove from list
                print_dbg(util::py::f("  Track {} with distance of {} is rejected with chi2 {:.3f}.", track_id, track_dist, track_chi2));
            }
        }

        return tracks_found_temp.size();
    }

    int VertexFinder::FindAll(std::vector<Track *> allTracks)
    {
        this->tracks_all = allTracks;

        if (tracks_all.size() < config["vertex_cut_VertexNHitsMin"])
            return 0;

        // Sort hits into groups
        print_dbg(util::py::f("\nVertex Finder: Event contains {} tracks", tracks_all.size()));

        // Seeding
        MakeSeeds();

        // Loop though required number of hits
        this->tracks_found_temp.clear();
        int ntracks_min = config["vertex_cut_VertexNHitsMin"];
        while (this->seeds_unused.size() > 0)
        {
            // Round 1: find hits based on seed
            int ntracks_found = FindOnceLS(seeds_unused.back());

            // Cut on number of tracks
            if (ntracks_found < ntracks_min)
            {
                print_dbg(util::py::f("-x Vertex rejected, only {} tracks", ntracks_found));
                seeds_unused.pop_back();
                continue;
            }

            // Print out the index of identified tracks
            std::string tracks_found_ids = "";
            for (auto track : tracks_found_temp)
            {
                tracks_found_ids += std::to_string(track->id) + " ";
            }
            print_dbg(util::py::f("-> Found vertex with {} tracks:", tracks_found_temp.size()), tracks_found_ids);

            // Round 2: Drop outliers
            // Discard this track if the seed contributes too much to chi2

            // Round 3: Run fit again using least square minimizer
            auto vertex_fitter = Kalman::LSVertex4DFitter(3, config["vertex_fit_MultipleScattering"]>0);
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
            std::unique_ptr<Tracker::Vertex> vertex_uniqueptr(vertex_found);
            this->vertices_found.push_back(std::move(vertex_uniqueptr));
            print_dbg(util::py::f("-> Vertex #{} added: chi2 {:.3f}", this->vertices_found.size(), this->vertices_found.back()->chi2));

            // Finally, remove other seeds that have hits in this vertex
            RemoveUsed();
            print_dbg("  Remaining seeds:", seeds_unused.size());
        }

        return vertices_found.size();
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

        info_nvertices = vertices_found.size();
        info_ntracks = 0;
        for (auto &vertex : vertices_found)
            info_ntracks += vertex->track_ids.size();

        return summary;
    }

} // namespace Tracker
