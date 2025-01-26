#include "track_finder.hh"
#include "kalman_model.hh"

namespace Tracker
{

    TrackFinder::TrackFinder(std::vector<DigiHit *> allHits,
                             bool debug, bool debug_kalman) : DEBUG(debug), DEBUG_KALMAN(debug_kalman), hits_all(allHits)
    {
        this->hits_found_temp = std::vector<DigiHit *>();
        this->tracks_found = TrackList();

        // Set the default parameters
        config["track_cut_SeedRange"] = 1;                                // [ns] Proper time of the seed pair
        config["track_cut_HitProjectionSigma"] = 10;                      // Range to look for the next hit, in unit of sigmas
        config["track_cut_HitAddChi2"] = 15;                              // Maximum chi2 for a hit to be added
        config["track_cut_HitDropChi2"] = 15;                             // Maximum chi2 for a hit to be kept during second rtrack_oun to -1 to turn off
        config["track_cut_TrackChi2Reduced"] = 3;                         // Cut on the chi2/ndof of the track
        config["track_cut_TrackChi2Prob"] = 0.9;                          // Cut on the  chi2 probability of the track
        config["track_cut_TrackNHitsMin"] = 4;                            // Cut on minimum number of hits
        config["track_cut_TrackSpeedLow"] = 25;                           // Cut on track speed [cm/ns]
        config["track_cut_TrackSpeedHigh"] = 35;                          // Cut on track speed [cm/ns]
        config["track_fit_MultipleScattering"] = 0;                       // Account for multiple scattering in fit. 0: OFF
        config["track_multiple_scattering_p"] = 500;                      // Particle momentum asumed for multiple scattering
        config["track_multiple_scattering_length"] = 0.06823501107481977; // Material thickness in the unit of scattering length. Calculated for 1 cm plastic scintillator and 0.6 mm aluminum.
    }

    void TrackFinder::Clear()
    {
    }

    void TrackFinder::MakeSeeds()
    {
        for (size_t i = 0; i < hits_all.size(); i++)
        {
            for (size_t j = i + 1; j < hits_all.size(); j++)
            {
                auto new_seed = new TrackSeed(hits_all[i], hits_all[j]);
                if (new_seed->score < config["track_cut_SeedRange"])
                    this->seeds_unused.push_back(new_seed);
                else
                    delete new_seed;
            }
        }
        // Sort seeds by the total distance from small to large. Larger distance should have higher priority.
        // Larger ones are at the end because pop_back is easier for vector
        std::sort(seeds_unused.begin(), seeds_unused.end(),
                  [](TrackSeed *a, TrackSeed *b) -> bool
                  { return a->dstep < b->dstep; });

        int iprint = 0;
        print_dbg(util::py::f("Making seeds, total {} seeds. Printing best 10. Smaller score is better.", seeds_unused.size()));
        if (DEBUG)
        {
            for (int i = seeds_unused.size() - 1; i >= 0; --i)
            {
                iprint += 1;
                auto seed = seeds_unused[i];
                print_dbg(util::py::f("  Seed [{}, {}]: score {:.4f}, dr {:.3f}, dt {:.3f}",
                                      seed->hits.first->id, seed->hits.second->id, seed->score, seed->dr, seed->dt));
                if (iprint > 10)
                    break;
            }
        }
    }

    void TrackFinder::GroupHits()
    {
        for (const auto &hit : this->hits_all)
        {
            if (hits_grouped.count(hit->group) == 0)
                hits_grouped[hit->group] = std::vector<DigiHit *>();
            hits_grouped[hit->group].push_back(hit);
        }

        // Write down all available groups
        this->hits_groups = std::vector<int>();
        for (const auto &pair : hits_grouped)
            this->hits_groups.push_back(pair.first);
    }

    int TrackFinder::FindOnce(TrackSeed *seed)
    {
        print_dbg(util::py::f("\nUsing Seed [{}, {}]", seed->hits.first->id, seed->hits.second->id));

        hits_found_temp.clear();
        auto finder = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, false);

        // Initialize the finder with seeds
        bool use_first_hit = true;
        hits_found_temp.push_back(seed->hits.first);
        hits_found_temp.push_back(seed->hits.second);
        finder.init_state(*seed->hits.first, *seed->hits.second, use_first_hit);

        // Loop all group except the one with the first hit
        for (const auto &pair : hits_grouped)
        {
            // Skip the group with hit in seed
            if (pair.first == seed->hits.first->group || pair.first == seed->hits.second->group)
                continue;
            auto hits_thisgroup = pair.second;
            print_dbg("  Searching in group", pair.first);

            // Update the internal matrices of finder with the position of this group
            finder.new_step(*hits_thisgroup[0]);

            // Calculate the chi2 contribution of each hit in this group
            // Only keep the minimum one
            int chi2_min_ind = -1;
            float chi2_min_val = std::numeric_limits<float>::infinity();

            for (size_t i = 0; i < hits_thisgroup.size(); i++)
            {
                auto hit = hits_thisgroup[i];
                float hit_chi2 = finder.try_measurement(*hit,
                                                        config["track_cut_HitProjectionSigma"],
                                                        config["track_cut_HitAddChi2"]);
                // print_dbg("    hit chi2:", hit_chi2);
                // find the smallest one
                if (hit_chi2 >= 0 && hit_chi2 < chi2_min_val)
                {
                    chi2_min_ind = i;
                    chi2_min_val = hit_chi2;
                }
            }
            print_dbg("    Minimum chi2 in this group is", chi2_min_val);

            // Add the hit with minimum chi2 and update the finder status
            if (chi2_min_ind != -1)
            {
                hits_found_temp.push_back(hits_thisgroup[chi2_min_ind]); // add to list
                finder.add_measurement(*hits_thisgroup[chi2_min_ind]);   // update finder
            }
        }

        return hits_found_temp.size();
    }

    int TrackFinder::FindAll()
    {
        // Seeding
        MakeSeeds();

        // Sort hits into groups
        GroupHits();
        print_dbg(util::py::f("Event contains hits in {} groups", hits_grouped.size()));

        if (hits_grouped.size() < config["track_cut_TrackNHitsMin"])
            return 0;

        // Loop though number of hits
        this->tracks_found.clear();
        int nhits_min = config["track_cut_TrackNHitsMin"];
        // for (int nhits_curr = hits_groups.size(); nhits_curr >= nhits_min; nhits_curr--)
        while (this->seeds_unused.size() > 0)
        {
            // Round 1: find hits based on seed
            int nhits_found = FindOnce(seeds_unused.back());
            seeds_unused.pop_back();
            // Cuts
            if (nhits_found < nhits_min)
            {
                print_dbg(util::py::f("-x Track rejected, not enough hits"));
                continue;
            }

            // Sort his by descending time.
            std::sort(hits_found_temp.begin(), hits_found_temp.end(),
                      [](DigiHit *a, DigiHit *b) -> bool
                      { return a->t() > b->t(); });

            // Debug print
            std::string hits_found_ids = "";
            for (auto hit: hits_found_temp)
                hits_found_ids += std::to_string(hit->id) + " ";
            print_dbg(util::py::f("-> Found track with {} hits:", nhits_found), hits_found_ids);


            // Round 3: run filter again backwards
            auto track_model = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, DEBUG_KALMAN);
            auto track_found = track_model.run_filter(hits_found_temp);
            this->tracks_found.push_back(std::move(track_found));
            print_dbg(util::py::f("-> Track #{} added: chi2 {:.3f}", this->tracks_found.size(), this->tracks_found.back()->chi2));

            // Finally, remove other seeds that have hits in this track
            RemoveUsed();
            print_dbg("  Remaining seeds:", seeds_unused.size());
        }

        return -1;
    }

    int TrackFinder::RemoveUsed()
    {
        // removing used seeds
        std::vector<TrackSeed *> seeds_unused_temp;
        std::vector<TrackSeed *> seeds_used_temp;
        std::string seeds_removed_list;
        for (size_t i = 0; i < seeds_unused.size(); i++)
        {
            bool veto_this = false;
            auto seed = seeds_unused[i];
            for (auto hit : hits_found_temp)
            {
                if (seed->hits.first->id == hit->id ||
                    seed->hits.second->id == hit->id)
                {
                    veto_this = true;
                    break;
                }
            }
            if (veto_this)
            {
                seeds_removed_list += "[" + std::to_string(seed->hits.first->id) + "," + std::to_string(seed->hits.second->id) + "] ";
                delete seed;
            }
            else
                seeds_unused_temp.push_back(seed);
        }
        print_dbg(util::py::f("  Seed removed: {}", seeds_removed_list));
        seeds_unused = seeds_unused_temp;

        // removing used hits
        for (auto hit : hits_found_temp)
        {
            for (int i = hits_grouped[hit->group].size() - 1; i >= 0; i--)
            {
                if (hits_grouped[hit->group][i]->id == hit->id)
                {
                    hits_grouped[hit->group].erase(hits_grouped[hit->group].begin() + i);
                }

            }
            if ( hits_grouped[hit->group].size()==0)
                hits_grouped.erase(hit->group);
        }

        return 0;
    }
} // namespace Tracker
