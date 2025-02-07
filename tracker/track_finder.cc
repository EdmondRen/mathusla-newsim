#include "track_finder.hh"

// ROOT libs
#include <TMath.h>
#include <Math/ProbFunc.h>

namespace Tracker
{

    TrackFinder::TrackFinder(
        bool debug, bool debug_kalman) : DEBUG(debug), DEBUG_KALMAN(debug_kalman)
    {
        this->hits_found_temp = std::vector<DigiHit *>();
        this->tracks_found = TrackList();

        // Set the default parameters
        config["track_cut_SeedRange"] = 1;                                // [ns] Proper time of the seed pair
        config["track_cut_HitProjectionSigma"] = 10;                      // Range to look for the next hit, in unit of sigmas
        config["track_cut_HitAddChi2First"] = 15;                         // Maximum chi2 for a hit to be added
        config["track_cut_HitAddChi2Other"] = 7;                          // Maximum chi2 for a hit to be added
        config["track_cut_HitDropChi2Prob"] = 0.9;                        // chi2 prob to drop a hit during smoothing
        config["track_cut_TrackChi2Reduced"] = 3;                         // Cut on the chi2/ndof of the track
        config["track_cut_TrackChi2Prob"] = 0.95;                         // Cut on the  chi2 probability of the track
        config["track_cut_TrackNHitsMin"] = 4;                            // Cut on minimum number of hits
        config["track_cut_TrackSpeedLow"] = 25;                           // Cut on track speed [cm/ns]
        config["track_cut_TrackSpeedHigh"] = 35;                          // Cut on track speed [cm/ns]
        config["track_fit_MultipleScattering"] = 1;                       // Account for multiple scattering in fit. 0: OFF
        config["track_multiple_scattering_p"] = 500;                      // Particle momentum asumed for multiple scattering
        config["track_multiple_scattering_length"] = 0.06823501107481977; // Material thickness in the unit of scattering length. Calculated for 1 cm plastic scintillator and 0.6 mm aluminum.
    }

    void TrackFinder::Config(std::map<std::string, double> &config_ext)
    {
        for (auto &pair : config)
        {
            if (config_ext.count(pair.first) > 0)
                config[pair.first] = config_ext[pair.first];
        }
    }

    void TrackFinder::Clear()
    {
    }

    void TrackFinder::MakeSeeds()
    {
        seeds_unused.clear();

        for (size_t i = 0; i < hits_all.size(); i++)
        {
            for (size_t j = i + 1; j < hits_all.size(); j++)
            {
                // Skip hits in the same or adjacent layer
                int dlayer = std::abs(hits_all[i]->layer - hits_all[j]->layer);
                if (dlayer == 0 || dlayer == 2)
                    continue;

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
                  {
                      if (a->dstep != b->dstep)
                      {
                          return a->dstep < b->dstep;
                      }
                      else
                          return a->dr < b->dr;
                  });

        int iprint = 0;
        print_dbg(util::py::f("Making seeds, total {} seeds. Printing best 30. Smaller score is better.", seeds_unused.size()));
        if (DEBUG)
        {
            for (int i = seeds_unused.size() - 1; i >= 0; --i)
            {
                iprint += 1;
                auto seed = seeds_unused[i];
                print_dbg(util::py::f("  Seed [{:03d}, {:03d}]: score {:.4f}, d_step: {:.0f}, dr {:.3f}, dt {:.3f}",
                                      seed->hits.first->id, seed->hits.second->id, seed->score, seed->dstep, seed->dr, seed->dt));
                if (iprint > 30)
                    break;
            }
        }
    }

    void TrackFinder::SortSeedsOccurance()
    {
        // std::unordered_map<int, int> hit_counts; // a map of {hit_id, hit_occurance}
        // for (auto &seed : seeds_unused)
        // {
        //     hit_counts[seed->hits.first->id] += 1;
        //     hit_counts[seed->hits.second->id] += 1;
        // }

        // for (auto &seed : seeds_unused)
        // {
        //     seed->nhit_occur += hit_counts[seed->hits.first->id];
        //     seed->nhit_occur += hit_counts[seed->hits.second->id];
        // }

        // std::sort(seeds_unused.begin(), seeds_unused.end(),
        //           [](TrackSeed *a, TrackSeed *b) -> bool
        //           {
        //               if (a->nhit_occur != b->nhit_occur)
        //               {
        //                   return a->nhit_occur > b->nhit_occur;
        //               }
        //               else if (a->dstep != b->dstep)
        //               {
        //                   return a->dstep < b->dstep;
        //               }
        //               else
        //                   return a->dr < b->dr;
        //           });

        // Sort seeds based on a commom direction
        // Vector3d direction = {0,0,0};
        // int counter = 0;
        // for (auto &seed : seeds_unused)
        // {
        //     direction += seed->dvec / seed->dr;
        //     counter += 1;
        // }
        // direction = direction/direction.norm();
        // print("Averaged direction", direction.transpose());
    }

    void TrackFinder::GroupHits()
    {
        hits_grouped.clear();

        for (const auto &hit : this->hits_all)
        {
            // if (hits_grouped.count(hit->layer) == 0)
            //     hits_grouped[hit->layer] = std::unordered_map<int, DigiHit *>();
            hits_grouped[hit->layer][hit->id] = hit;
        }

        // Write down all available groups
        this->hits_groups = std::vector<int>();
        for (const auto &pair : hits_grouped)
            this->hits_groups.push_back(pair.first);
    }

    int TrackFinder::FindOnce(TrackSeed *seed)
    {
        print_dbg(util::py::f("Using Seed [{},{}]", seed->hits.first->id, seed->hits.second->id));

        hits_found_temp.clear();
        auto finder = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, false);

        // Initialize the finder with seeds
        bool use_first_hit = true;
        hits_found_temp.push_back(seed->hits.first);
        hits_found_temp.push_back(seed->hits.second);
        finder.init_state(*seed->hits.first, *seed->hits.second, use_first_hit);

        // Loop all layers except the one with the first hit
        for (const auto &pair : hits_grouped)
        {
            // pair: {layer, hits}
            // Skip the layer with hit in seed
            if (pair.first == seed->hits.first->layer || pair.first == seed->hits.second->layer)
                continue;
            auto hits_thisgroup = pair.second;

            // Update the internal matrices of finder with the position of this layer
            finder.new_step(*hits_thisgroup.begin()->second);

            // Calculate the chi2 contribution of each hit in this layer
            // Only keep the minimum one
            int chi2_min_ind = -1;
            float chi2_min_val = std::numeric_limits<float>::infinity();

            for (const auto &hit_and_id : hits_thisgroup)
            {
                auto hit_id = hit_and_id.first;
                auto hit = hit_and_id.second;
                float hit_chi2;
                if (hits_found_temp.size() == 2)
                    hit_chi2 = finder.try_measurement(*hit,
                                                      config["track_cut_HitProjectionSigma"],
                                                      config["track_cut_HitAddChi2First"]);
                else
                    hit_chi2 = finder.try_measurement(*hit,
                                                      config["track_cut_HitProjectionSigma"],
                                                      config["track_cut_HitAddChi2Other"]);
                // print_dbg("    hit chi2:", hit_chi2);
                // find the smallest one
                if (hit_chi2 >= 0 && hit_chi2 < chi2_min_val)
                {
                    chi2_min_ind = hit_id;
                    chi2_min_val = hit_chi2;
                }
            }
            print_dbg(util::py::f("  Searching layer {}. Minimum chi2 is {:.3f}, from hit {}.", pair.first, chi2_min_val, chi2_min_ind));

            // Add the hit with minimum chi2 and update the finder status
            if (chi2_min_ind != -1)
            {
                hits_found_temp.push_back(hits_thisgroup[chi2_min_ind]); // add to list
                finder.add_measurement(*hits_thisgroup[chi2_min_ind]);   // update finder
            }
        }

        hits_found_temp_chi2 = finder.get_current_chi2();

        return hits_found_temp.size();
    }

    int TrackFinder::FindAll(std::vector<DigiHit *> allHits)
    {
        this->hits_all = allHits;

        // Sort hits into groups
        GroupHits();
        print_dbg(util::py::f("\nTrack Finder: Event contains {} hits in {} groups", hits_all.size(), hits_grouped.size()));

        // Seeding
        MakeSeeds();

        if (hits_grouped.size() < config["track_cut_TrackNHitsMin"])
            return 0;

        // Loop though required number of hits
        this->tracks_found.clear();
        int nhits_min = config["track_cut_TrackNHitsMin"];
        for (int nhits_curr = hits_groups.size(); nhits_curr >= nhits_min; nhits_curr--)
        {
            print_dbg(util::py::f("\n***Searching for tracks with at least {} hits***", nhits_curr));

            //////////////
            if (nhits_curr == 4)
            {
                // SortSeedsOccurance();
                // break;
            }

            while (this->seeds_unused.size() > 0)
            {
                // Skip the seed that have been used and known to have fewer hits than required
                if (seeds_unused.back()->nhits_found > 0 && seeds_unused.back()->nhits_found < nhits_curr)
                {
                    seeds_unused_next.insert(seeds_unused_next.begin(), seeds_unused.back());
                    seeds_unused.pop_back();
                    continue;
                }

                // **************************************************************************
                // Round 1: find hits based on seed
                int nhits_found = FindOnce(seeds_unused.back());

                // Cut on number of hits
                if (nhits_found < nhits_min)
                {
                    print_dbg(util::py::f("-x Track rejected, only {} hits", nhits_found));
                    seeds_unused.pop_back();
                    continue;
                }

                // Sort hits by descending time.
                std::sort(hits_found_temp.begin(), hits_found_temp.end(),
                          [](DigiHit *a, DigiHit *b) -> bool
                          { return a->t() > b->t(); });
                std::string hits_found_ids = "";
                for (auto hit : hits_found_temp)
                {
                    hits_found_ids += std::to_string(hit->id) + " ";
                }

                // Skip to next round if nhits is fewer than asked
                if (nhits_found < nhits_curr)
                {
                    seeds_unused.back()->nhits_found = nhits_found;
                    seeds_unused.back()->chi2prob_found = ROOT::Math::chisquared_cdf(hits_found_temp_chi2, 3 * nhits_found - 6);
                    if (seeds_unused.back()->chi2prob_found > config["track_cut_TrackChi2Prob"])
                    {
                        print_dbg("-x Track rejected, only {} hits and chi2 is larger than specified");
                    }
                    else
                    {
                        seeds_unused_next.insert(seeds_unused_next.begin(), seeds_unused.back());
                        print_dbg(util::py::f("-x Track rejected, only {} hits. Keep the seed for next iteration. hit inds:", nhits_found), hits_found_ids);
                    }

                    seeds_unused.pop_back();
                    continue;
                }

                print_dbg(util::py::f("-> Found track with {} hits:", nhits_found), hits_found_ids);

                // **************************************************************************
                // Round 2: drop hits with excessive chi2
                auto track_model = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, DEBUG_KALMAN); // DEBUG_KALMAN

                if (config["track_cut_HitDropChi2Prob"] > 0)
                {
                    auto track_model_drop = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, DEBUG_KALMAN); // DEBUG_KALMAN
                    track_model_drop.run_filter_smooth(hits_found_temp, config["track_cut_HitDropChi2Prob"]);
                    auto dropped_inds = track_model_drop.get_dropped_inds();
                    if (dropped_inds.size() > 0)
                    {
                        for (auto ind : dropped_inds)
                        {
                            print_dbg(util::py::f("-x Hit {} dropped during smoothing due to large chi2", hits_found_temp[ind]->id));
                            hits_found_temp.erase(hits_found_temp.begin() + ind);
                        }
                    }
                    if (static_cast<int>(hits_found_temp.size()) < nhits_min)
                    {
                        print_dbg(util::py::f("-x Track rejected, only {} hits after dropping", hits_found_temp.size()));
                        seeds_unused.pop_back();
                        continue;
                    }
                }

                // std::vector<double> hits_chi2;
                // for (auto hit : hits_found_temp)
                // {
                //     auto pos_and_cov = track_temp->get_same_invar_pos_and_cov(hit->vec4, -1, config["track_fit_MultipleScattering"] > 0);
                //     // Residual
                //     Vector4d residual4 = pos_and_cov.first - hit->vec4;
                //     auto residual = Helper::removeElement(residual4, track_temp->iv_index);
                //     // Cov
                //     auto cov_mod = pos_and_cov.second;
                //     cov_mod.diagonal() += hit->vec3_err.array().square().matrix();
                //     // chi2
                //     double chi2_temp = residual.transpose() * cov_mod.inverse() * residual;
                //     hits_chi2.push_back(chi2_temp);
                // }

                // track_model.new_step(*seeds_unused.back()->hits.first);
                // auto chi2_seed1 = track_model.try_measurement(*seeds_unused.back()->hits.first, 5, 1000);

                // track_model.new_step(*seeds_unused.back()->hits.second);
                // auto chi2_seed2 = track_model.try_measurement(*seeds_unused.back()->hits.second, 5, 1000);
                // print_dbg("  chi2 of seed hits:", chi2_seed1, chi2_seed2);
                // if (chi2_seed1 > config["track_cut_HitDropChi2"] || chi2_seed2 > config["track_cut_HitDropChi2"])
                // {
                //     print_dbg("-x hit from seed have large chi2. Probablity wrong combination. Track rejected.");
                //     seeds_unused.pop_back();
                //     continue;
                // }

                // **************************************************************************
                // Round 3: Run filter again backwards
                // Discard this track if the seed contributes too much to chi2
                // auto track_model = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, DEBUG_KALMAN); // DEBUG_KALMAN
                auto track_found = track_model.run_filter(hits_found_temp);

                // Cut on chi2
                int ndof = hits_found_temp.size() * 3 - 6;
                float chi2 = track_found->chi2;
                float chi2_prob = ROOT::Math::chisquared_cdf(chi2, ndof);
                // bool passed = (ndof == 6 && (chi2 / ndof) < config["track_cut_TrackChi2Reduced"]) ||
                //               (ndof > 6 && chi2_prob < config["track_cut_TrackChi2Prob"]);
                bool passed = chi2_prob < config["track_cut_TrackChi2Prob"];
                if (!passed)
                {
                    print_dbg(util::py::f("-x Track rejected, chi2/ndof ({}/{}) exceeds limit, prob = {}", chi2, ndof, 1 - chi2_prob));
                    seeds_unused.pop_back();
                    continue;
                }

                // Setting the metadata of this track
                track_found->id = this->tracks_found.size();
                track_found->hit_ids = std::vector<int>();
                track_found->convert_to_time(); // Add a set of track parameter using time as independent variable.

                hits_found_ids = "";
                for (auto hit : hits_found_temp)
                {
                    hits_found_ids += std::to_string(hit->id) + " ";
                    track_found->hit_ids.push_back(hit->id);
                }
                print_dbg(util::py::f("-> Found track with {} hits:", nhits_found), hits_found_ids);

                // Add track to the list
                this->tracks_found.push_back(std::move(track_found));
                print_dbg(util::py::f("-> Track #{} added: chi2 {:.3f}", this->tracks_found.size(), this->tracks_found.back()->chi2));

                // Finally, remove other seeds that have hits in this track
                RemoveUsed();
                print_dbg("  Remaining seeds:", seeds_unused.size());
            }

            // Sort the remaining seed, and reset the unused seed pointer
            std::sort(this->seeds_unused_next.begin(), this->seeds_unused_next.end(),
                      [](auto *a, auto *b)
                      {
                          if (a->nhits_found != b->nhits_found)
                          {
                              return a->nhits_found < b->nhits_found;
                          }
                          else
                              return a->chi2prob_found > b->chi2prob_found;
                      });

            this->seeds_unused = std::move(this->seeds_unused_next);

            // // Let's have a look at the situation of 4 hits
            // if (nhits_curr == 5)
            // {
            //     print("Remaining seeds");
            //     int iprint = 0;
            //     for (int i = seeds_unused.size() - 1; i >= 0; --i)
            //     {
            //         iprint += 1;
            //         auto seed = seeds_unused[i];
            //         print_dbg(util::py::f("  Seed [{:03d}, {:03d}]: chi2 {:.4f}, d_step: {:.0f}, dr {:.3f}, dt {:.3f}",
            //                               seed->hits.first->id, seed->hits.second->id, seed->chi2prob_found, seed->dstep, seed->dr, seed->dt));
            //         // if (iprint > 30)
            //         //     break;
            //     }

            //     // Loop all layers except the one with the first hit
            //     for (const auto &pair : hits_grouped)
            //     {
            //         // pair: {layer, hits}

            //         auto hits_thisgroup = pair.second;
            //             print("Remaing hits in layer", pair.first);

            //         for (const auto &hit_and_id : hits_thisgroup)
            //         {
            //             auto hit_id = hit_and_id.first;
            //             auto hit = hit_and_id.second;
            //             print("Hit ", hit_id);
            //         }
            //     }
            // }
        }

        // for (auto &seed : seeds_unused)
        // {

        //     // Round 1: find hits based on seed
        //     int nhits_found = FindOnce(seed);

        //     // Cut on number of hits
        //     if (nhits_found < nhits_min)
        //     {
        //         print_dbg(util::py::f("-x Track rejected, only {} hits", nhits_found));
        //         seeds_unused.pop_back();

        //         continue;
        //     }

        //     // Sort hits by descending time.
        //     std::sort(hits_found_temp.begin(), hits_found_temp.end(),
        //               [](DigiHit *a, DigiHit *b) -> bool
        //               { return a->t() > b->t(); });
        //     std::string hits_found_ids = "";
        //     for (auto hit : hits_found_temp)
        //     {
        //         hits_found_ids += std::to_string(hit->id) + " ";
        //     }

        //     // Skip to next round if nhits is fewer than asked
        //     if (nhits_found < 4)
        //     {
        //         print_dbg(util::py::f("-x Track rejected, only {} hits. Keep the seed for next iteration. hit inds:", nhits_found), hits_found_ids);
        //         continue;
        //     }
        //     print_dbg(util::py::f("-> Found track with {} hits:", nhits_found), hits_found_ids,"chi2:", this->hits_found_temp_chi2);
        // }
        // if (seeds_unused.size() > 2)
        // {
        //     for (int i = 0; i < seeds_unused.size() - 1; i++)
        //     {
        //         for (int j = i + 1; j < seeds_unused.size(); j++)
        //         {
        //             auto seed1 = seeds_unused[i];
        //             auto seed2 = seeds_unused[j];
        //             if (seed1->hits.first->layer == seed2->hits.first->layer ||
        //                 seed1->hits.first->layer == seed2->hits.second->layer ||
        //                 seed1->hits.second->layer == seed2->hits.first->layer ||
        //                 seed1->hits.second->layer == seed2->hits.second->layer)
        //                 continue;
        //             std::vector<DigiHit *> hits_try = {seed1->hits.first, seed1->hits.second, seed2->hits.first, seed2->hits.second};
        //             std::sort(hits_try.begin(), hits_try.end(),
        //                       [](DigiHit *a, DigiHit *b) -> bool
        //                       { return a->t() > b->t(); });
        //             auto track_model = Kalman::KalmanTrack4D(config["track_fit_MultipleScattering"], 4, 6, DEBUG_KALMAN); // DEBUG_KALMAN
        //             auto track = track_model.run_filter(hits_try);
        //             if (track->chi2 < 18)
        //                 print_dbg(util::py::f("-> Found track with chi2 {} : hits:", track->chi2), seed1->hits.first->id, seed1->hits.second->id, seed2->hits.first->id, seed2->hits.second->id);
        //         }
        //     }
        // }

        if (DEBUG)
        {
            print("\n\n------------------------------------------------");
            auto info = this->Summary();
            print(info);

            for (const auto &track : tracks_found)
            {
                print(util::py::f(" * Track #{}, nhits {}, chi2 {}, params = ", track->id, track->hit_ids.size(), track->chi2), track->params.transpose());
            }
            print("------------------------------------------------\n\n");
        }

        return 0;
    }

    int TrackFinder::RemoveUsed()
    {
        // removing used seeds
        std::vector<TrackSeed *> seeds_unused_temp;
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

        // removing used seeds
        std::vector<TrackSeed *> seeds_unused_next_temp;
        seeds_removed_list = "";
        for (size_t i = 0; i < seeds_unused_next.size(); i++)
        {
            bool veto_this = false;
            auto seed = seeds_unused_next[i];
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
                seeds_unused_next_temp.push_back(seed);
        }
        print_dbg(util::py::f("  Seed removed: {}", seeds_removed_list));
        seeds_unused_next = seeds_unused_next_temp;

        // removing used hits
        for (auto hit : hits_found_temp)
        {
            hits_grouped[hit->layer].erase(hit->id);
            if (hits_grouped[hit->layer].size() == 0)
                hits_grouped.erase(hit->layer);
        }

        return 0;
    }

    std::string TrackFinder::Summary()
    {

        summary = "";
        summary += util::py::f("Event track-finding summary:\n");
        summary += util::py::f("  Total hits: {}\n", hits_all.size());
        summary += util::py::f("  Total tracks: {}\n", tracks_found.size());

        int nhits_min = config["track_cut_TrackNHitsMin"];
        for (int nhits_curr = hits_groups.size(); nhits_curr >= nhits_min; nhits_curr--)
        {
            int ntrack = 0;
            for (const auto &track : tracks_found)
            {
                if (static_cast<int>(track->hit_ids.size()) == nhits_curr)
                    ntrack += 1;
            }
            summary += util::py::f("   # tracks with {} hits: {}\n", nhits_curr, ntrack);
        }

        info_ntracks = tracks_found.size();
        info_nhits = 0;
        for (const auto &track : tracks_found)
            info_nhits += track->hit_ids.size();

        return summary;
    }

    // ******************************************************************************************************************

    VertexTrackFinder::VertexTrackFinder(
        bool debug, bool debug_kalman) : DEBUG(debug), DEBUG_KALMAN(debug_kalman)
    {
        this->hits_found_temp = std::vector<DigiHit *>();
        this->tracks_found = TrackList();

        // Set the default parameters
        config["tracklet_cut_SeedRange"] = 1;           // [ns] Proper time of the seed pair
        config["tracklet_cut_HitProjectionSigma"] = 10; // Range to look for the next hit, in unit of sigmas
        config["tracklet_cut_HitAddChi2First"] = 15;    // Maximum chi2 for a hit to be added
        config["tracklet_cut_HitAddChi2Other"] = 7;     // Maximum chi2 for a hit to be added
        config["tracklet_cut_TrackChi2Prob"] = 0.95;    // Cut on the  chi2 probability of the track
        config["tracklet_cut_TrackNHitsMin"] = 2;       // Cut on minimum number of hits
        config["tracklet_fit_MultipleScattering"] = 1;  // Account for multiple scattering in fit. 0: OFF

        finder = new Kalman::KalmanTrack4D(config["tracklet_fit_MultipleScattering"], 4, 6, false);
    }

    void VertexTrackFinder::Config(std::map<std::string, double> &config_ext)
    {
        for (auto &pair : config)
        {
            if (config_ext.count(pair.first) > 0)
                config[pair.first] = config_ext[pair.first];
        }
    }

    void VertexTrackFinder::Clear()
    {
    }

    void VertexTrackFinder::MakeSeeds()
    {
        seeds_unused.clear();

        for (auto igroup : this->hits_groups)
        {
            for (auto hit_pair : hits_grouped[igroup])
            {
                auto new_seed = new VertexTrackSeed(this->vertex->params, hit_pair.second);
                new_seed->layer = hit_pair.second->layer;
                if (new_seed->score < config["tracklet_cut_SeedRange"])
                    this->seeds_unused.push_back(new_seed);
                else
                    delete new_seed;
            }
        }
        // Sort seeds by the total distance from small to large. Larger distance should have higher priority.
        // Larger ones are at the end because pop_back is easier for vector
        std::sort(seeds_unused.begin(), seeds_unused.end(),
                  [](VertexTrackSeed *a, VertexTrackSeed *b) -> bool
                  {
                      if (a->layer != b->layer)
                          return a->layer > b->layer;
                      else
                          return a->dr < b->dr;
                  });

        int iprint = 0;
        print_dbg(util::py::f("Making seeds, total {} seeds. Printing best 30. Smaller score is better.", seeds_unused.size()));
        if (DEBUG)
        {
            for (int i = seeds_unused.size() - 1; i >= 0; --i)
            {
                iprint += 1;
                auto seed = seeds_unused[i];
                print_dbg(util::py::f("  Seed [{:03d}]: score {:.4f}, dr {:.3f}, dt {:.3f}",
                                      seed->hits.first, seed->score, seed->dr, seed->dt));
                if (iprint > 30)
                    break;
            }
        }
    }

    void VertexTrackFinder::GroupHits()
    {
        hits_grouped.clear();

        // Get list of hits already used by the vertex
        std::unordered_map<int, int> hits_vertex_inds;
        for (auto i : vertex->hit_ids)
        {
            hits_vertex_inds[i] = 0;
        }

        for (const auto &hit : this->hits_all)
        {
            // Use only hits that are not part of the vertex
            if (hits_vertex_inds.count(hit->id) > 0)
                hits_grouped[hit->layer][hit->id] = hit;
        }

        // Write down all available groups
        this->hits_groups = std::vector<int>();
        for (const auto &pair : hits_grouped)
            this->hits_groups.push_back(pair.first);
        std::sort(this->hits_groups.begin(), this->hits_groups.end());
    }

    int VertexTrackFinder::FindOnce(VertexTrackSeed *seed)
    {
        print_dbg(util::py::f("Using Seed [{}]", seed->hits.first));

        hits_found_temp.clear();

        // Initialize the finder with seeds
        hits_found_temp.push_back(seed->hits.second);
        finder->init_state(this->vertex->params, *seed->hits.second, this->vertex->cov);

        // Loop all layers except the one with the first hit
        for (const auto &pair : hits_grouped)
        {
            // pair: {layer, hits}
            // Skip the layer with hit in seed
            if (pair.first == seed->hits.second->layer)
                continue;
            auto hits_thisgroup = pair.second;

            // Update the internal matrices of finder with the position of this layer
            finder->new_step(*hits_thisgroup.begin()->second);

            // Calculate the chi2 contribution of each hit in this layer
            // Only keep the minimum one
            int chi2_min_ind = -1;
            float chi2_min_val = std::numeric_limits<float>::infinity();

            for (const auto &hit_and_id : hits_thisgroup)
            {
                auto hit_id = hit_and_id.first;
                auto hit = hit_and_id.second;
                float hit_chi2;
                if (hits_found_temp.size() == 2)
                    hit_chi2 = finder->try_measurement(*hit,
                                                       config["tracklet_cut_HitProjectionSigma"],
                                                       config["tracklet_cut_HitAddChi2First"]);
                else
                    hit_chi2 = finder->try_measurement(*hit,
                                                       config["tracklet_cut_HitProjectionSigma"],
                                                       config["tracklet_cut_HitAddChi2Other"]);
                // print_dbg("    hit chi2:", hit_chi2);
                // find the smallest one
                if (hit_chi2 >= 0 && hit_chi2 < chi2_min_val)
                {
                    chi2_min_ind = hit_id;
                    chi2_min_val = hit_chi2;
                }
            }
            print_dbg(util::py::f("  Searching layer {}. Minimum chi2 is {:.3f}, from hit {}.", pair.first, chi2_min_val, chi2_min_ind));

            // Add the hit with minimum chi2 and update the finder status
            if (chi2_min_ind != -1)
            {
                hits_found_temp.push_back(hits_thisgroup[chi2_min_ind]); // add to list
                finder->add_measurement(*hits_thisgroup[chi2_min_ind]);  // update finder
            }
        }

        hits_found_temp_chi2 = finder->get_current_chi2();

        return hits_found_temp.size();
    }

    int VertexTrackFinder::FindAll(std::vector<DigiHit *> allHits, Vertex *_vertex_)
    {
        this->hits_all = allHits;
        this->vertex = _vertex_;

        // Sort hits into groups
        GroupHits();
        print_dbg(util::py::f("\nTrack Finder: Event contains {} hits in {} groups", hits_all.size(), hits_grouped.size()));

        // Seeding
        MakeSeeds();

        if (hits_grouped.size() < config["tracklet_cut_TrackNHitsMin"])
            return 0;

        // Loop though required number of hits
        this->tracks_found.clear();
        int nhits_min = config["tracklet_cut_TrackNHitsMin"];
        for (int nhits_curr = hits_groups.size(); nhits_curr >= nhits_min; nhits_curr--)
        {
            print_dbg(util::py::f("\n***Searching for tracks with at least {} hits***", nhits_curr));

            while (this->seeds_unused.size() > 0)
            {
                // Skip the seed that have been used and known to have fewer hits than required
                if (seeds_unused.back()->nhits_found > 0 && seeds_unused.back()->nhits_found < nhits_curr)
                {
                    seeds_unused_next.insert(seeds_unused_next.begin(), seeds_unused.back());
                    seeds_unused.pop_back();
                    continue;
                }

                // **************************************************************************
                // Round 1: find hits based on seed
                int nhits_found = FindOnce(seeds_unused.back());

                // Cut on number of hits
                if (nhits_found < nhits_min)
                {
                    print_dbg(util::py::f("-x Track rejected, only {} hits", nhits_found));
                    seeds_unused.pop_back();
                    continue;
                }

                // Sort hits by descending time.
                std::sort(hits_found_temp.begin(), hits_found_temp.end(),
                          [](DigiHit *a, DigiHit *b) -> bool
                          { return a->t() > b->t(); });
                std::string hits_found_ids = "";
                for (auto hit : hits_found_temp)
                {
                    hits_found_ids += std::to_string(hit->id) + " ";
                }

                // Skip to next round if nhits is fewer than asked
                int ndof = nhits_found * 3 - 6;
                if (nhits_found < nhits_curr)
                {
                    seeds_unused.back()->nhits_found = nhits_found;
                    seeds_unused.back()->chi2prob_found = ROOT::Math::chisquared_cdf(hits_found_temp_chi2, ndof);
                    if (seeds_unused.back()->chi2prob_found > config["tracklet_cut_TrackChi2Prob"])
                    {
                        print_dbg("-x Track rejected, only {} hits and chi2 is larger than specified");
                    }
                    else
                    {
                        seeds_unused_next.insert(seeds_unused_next.begin(), seeds_unused.back());
                        print_dbg(util::py::f("-x Track rejected, only {} hits. Keep the seed for next iteration. hit inds:", nhits_found), hits_found_ids);
                    }

                    seeds_unused.pop_back();
                    continue;
                }

                print_dbg(util::py::f("-> Found track with {} hits:", nhits_found), hits_found_ids);

                // Cut on chi2
                float chi2 = hits_found_temp_chi2;
                float chi2_prob = ROOT::Math::chisquared_cdf(chi2, ndof);
                bool passed = chi2_prob < config["tracklet_cut_TrackChi2Prob"];
                if (!passed)
                {
                    print_dbg(util::py::f("-x Track rejected, chi2/ndof ({}/{}) exceeds limit, prob = {}", chi2, ndof, 1 - chi2_prob));
                    seeds_unused.pop_back();
                    continue;
                }

                // Make the track
                std::unique_ptr<Track> track_found = std::make_unique<Track>();
                track_found->params = finder->kf.GetState();
                track_found->cov = finder->kf.GetCov();
                track_found->chi2 = finder->kf.GetChi2();
                track_found->iv_index = hits_found_temp.back()->iv_index;
                track_found->iv_value = hits_found_temp.back()->get_step();
                track_found->iv_error = hits_found_temp.back()->get_steperr();

                // Setting the metadata of this track
                track_found->id = this->tracks_found.size();
                track_found->hit_ids = std::vector<int>();

                hits_found_ids = "";
                for (auto hit : hits_found_temp)
                {
                    hits_found_ids += std::to_string(hit->id) + " ";
                    track_found->hit_ids.push_back(hit->id);
                }
                print_dbg(util::py::f("-> Found track with {} hits:", nhits_found), hits_found_ids);

                // Add track to the list
                this->tracks_found.push_back(std::move(track_found));
                print_dbg(util::py::f("-> Track #{} added: chi2 {:.3f}", this->tracks_found.size(), this->tracks_found.back()->chi2));

                // Finally, remove other seeds that have hits in this track
                RemoveUsed();
                print_dbg("  Remaining seeds:", seeds_unused.size());
            }

            // Sort the remaining seed, and reset the unused seed pointer
            std::sort(this->seeds_unused_next.begin(), this->seeds_unused_next.end(),
                      [](auto *a, auto *b)
                      {
                          if (a->nhits_found != b->nhits_found)
                          {
                              return a->nhits_found < b->nhits_found;
                          }
                          else
                              return a->chi2prob_found > b->chi2prob_found;
                      });

            this->seeds_unused = std::move(this->seeds_unused_next);
        }

        if (DEBUG)
        {
            print("\n\n------------------------------------------------");
            auto info = this->Summary();
            print(info);

            for (const auto &track : tracks_found)
            {
                print(util::py::f(" * Track #{}, nhits {}, chi2 {}, params = ", track->id, track->hit_ids.size(), track->chi2), track->params.transpose());
            }
            print("------------------------------------------------\n\n");
        }

        return 0;
    }

    int VertexTrackFinder::RemoveUsed()
    {
        // removing used seeds
        std::vector<VertexTrackSeed *> seeds_unused_temp;
        std::string seeds_removed_list;
        for (size_t i = 0; i < seeds_unused.size(); i++)
        {
            bool veto_this = false;
            auto seed = seeds_unused[i];
            for (auto hit : hits_found_temp)
            {
                if (seed->hits.second->id == hit->id)
                {
                    veto_this = true;
                    break;
                }
            }
            if (veto_this)
            {
                seeds_removed_list += "[" + std::to_string(seed->hits.second->id) + "] ";
                delete seed;
            }
            else
                seeds_unused_temp.push_back(seed);
        }
        print_dbg(util::py::f("  Seed removed: {}", seeds_removed_list));
        seeds_unused = seeds_unused_temp;

        // removing used seeds
        std::vector<VertexTrackSeed *> seeds_unused_next_temp;
        seeds_removed_list = "";
        for (size_t i = 0; i < seeds_unused_next.size(); i++)
        {
            bool veto_this = false;
            auto seed = seeds_unused_next[i];
            for (auto hit : hits_found_temp)
            {
                if (seed->hits.second->id == hit->id)
                {
                    veto_this = true;
                    break;
                }
            }
            if (veto_this)
            {
                seeds_removed_list += "[" + std::to_string(seed->hits.second->id) + "] ";
                delete seed;
            }
            else
                seeds_unused_next_temp.push_back(seed);
        }
        print_dbg(util::py::f("  Seed removed: {}", seeds_removed_list));
        seeds_unused_next = seeds_unused_next_temp;

        // removing used hits
        for (auto hit : hits_found_temp)
        {
            hits_grouped[hit->layer].erase(hit->id);
            if (hits_grouped[hit->layer].size() == 0)
                hits_grouped.erase(hit->layer);
        }

        return 0;
    }

    std::string VertexTrackFinder::Summary()
    {

        summary = "";
        summary += util::py::f("Event track-finding summary:\n");
        summary += util::py::f("  Total hits: {}\n", hits_all.size());
        summary += util::py::f("  Total tracks: {}\n", tracks_found.size());

        int nhits_min = config["tracklet_cut_TrackNHitsMin"];
        for (int nhits_curr = hits_groups.size(); nhits_curr >= nhits_min; nhits_curr--)
        {
            int ntrack = 0;
            for (const auto &track : tracks_found)
            {
                if (static_cast<int>(track->hit_ids.size()) == nhits_curr)
                    ntrack += 1;
            }
            summary += util::py::f("   # tracks with {} hits: {}\n", nhits_curr, ntrack);
        }

        info_ntracks = tracks_found.size();
        info_nhits = 0;
        for (const auto &track : tracks_found)
            info_nhits += track->hit_ids.size();

        return summary;
    }

} // namespace Tracker
