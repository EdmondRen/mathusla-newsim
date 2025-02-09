#include <math.h>
#include "types_and_helper.hh"

namespace Tracker
{
    // -------------------------------------------------------------------------------------------
    // Class Track

    VectorXd Track::convert_to_time()
    {
        if (iv_index < 3)
        {
            params_time = VectorXd(6);
            cov_time = MatrixXd(6, 6);
            t0 = params(2);
            MatrixXdSp jac(6, 6);

            // Setting parameters
            for (int j = 0, k = 0; j < 3; ++j)
            {
                if (j == iv_index)
                {
                    params_time(j) = iv_value;
                    params_time(j + 3) = 1 / params(2 + 3);
                }
                else
                {
                    params_time(j) = params(k);
                    params_time(j + 3) = params(k + 3) / params(2 + 3);
                    k++;
                }
            }

            // Setting Jacobian
            this->t0_err = std::pow(cov(2, 2), 0.5);
            for (int j = 0; j < 3; ++j)
            {
                if (j == iv_index)
                {
                    jac.insert(j, j) = params_time(j + 3);
                    jac.insert(j + 3, j + 3) = -params_time(j + 3) * params_time(j + 3);
                    jac.insert(j, j + 3) = -params_time(j + 3) * params_time(j + 3) * t0_err;
                }
                else
                {
                    jac.insert(j, j) = 1;
                    jac.insert(j, iv_index) = params_time(j + 3);
                    jac.insert(j + 3, j + 3) = params_time(iv_index + 3);
                    jac.insert(j + 3, iv_index + 3) = -params_time(j + 3) * params_time(iv_index + 3);
                    jac.insert(j, j + 3) = params_time(iv_index + 3) * t0_err;
                    jac.insert(j, iv_index + 3) = -params_time(j + 3) * params_time(iv_index + 3) * t0_err;
                }
            }
            cov_time = jac * cov * jac.transpose();
            // std::cout << "Jac \n"
            //           << jac << std::endl;
            // std::cout << "cov original \n"
            //           << cov << std::endl;
            // std::cout << "cov time \n"
            //           << cov_time << std::endl;
        }
        else
        {
            params_time = params; // Do nothing if it is already using time as independent variable.
            cov_time = cov;
            t0 = iv_value;
            t0_err = iv_error;
        }

        return params_time;
    }

    std::pair<double, Vector4d> Track::get_closest_midpoint(const Track &track1, const Track &track2)
    {
        VectorXd s1 = track1.params_time;
        VectorXd s2 = track2.params_time;

        auto r1 = s1.segment(0, 3);
        auto r2 = s2.segment(0, 3);
        auto v1 = s1.segment(3, 3);
        auto v2 = s2.segment(3, 3);
        auto t1 = track1.t0;
        auto t2 = track2.t0;
        auto dr = r1 - r2;
        auto dv = v2 - v1;
        double time_midpoint = (dr.dot(dv) + (v2 * t2 - v1 * t1).dot(dv)) / std::pow(dv.norm(), 2);

        Vector3d pos1 = r1 + v1 * (time_midpoint - t1);
        Vector3d pos2 = r2 + v2 * (time_midpoint - t2);
        Vector3d pos_midpoint = 0.5 * (pos1 + pos2);
        double distance = (pos1 - pos2).norm();

        Vector4d midpoint4d(pos_midpoint(0), pos_midpoint(1), pos_midpoint(2), time_midpoint);

        return std::make_pair(distance, midpoint4d);
    }

    double Track::get_same_time_dist(Vector4d point) const
    {
        double dt = point(3) - this->t0;
        Vector3d pos_on_track = this->params_time.segment(0, 3) + this->params_time.segment(3, 3) * dt;
        Vector3d pos_residual = pos_on_track - point.segment(0, 3);
        double dist_point = pos_residual.norm();

        return dist_point;
    }

    std::pair<Vector4d, MatrixXd> Track::get_same_time_pos_and_cov(Vector4d point, double speed_constraint, bool multiple_scattering) const
    {
        (void)speed_constraint;
        (void)multiple_scattering;

        double dt = point(3) - this->t0;
        Vector3d pos_on_track = this->params_time.segment(0, 3) + this->params_time.segment(3, 3) * dt;
        Vector4d pos_full = Helper::insertVector(pos_on_track, 3, point(3));
        MatrixXd cov_residual = this->cov_time.topLeftCorner(3, 3) +
                                this->cov_time.topRightCorner(3, 3) * dt +
                                this->cov_time.bottomLeftCorner(3, 3) * dt +
                                this->cov_time.bottomRightCorner(3, 3) * dt * dt;

        return std::make_pair(pos_full, cov_residual);
    }

    std::pair<double, double> Track::get_same_time_dist_and_chi2(Vector4d point, double speed_constraint, bool multiple_scattering) const
    {
        (void)speed_constraint;

        auto point_and_cov = get_same_time_pos_and_cov(point, speed_constraint, multiple_scattering);
        Vector3d residual = point_and_cov.first.segment(0, 3) - point.segment(0, 3);
        double dist_point = residual.norm();
        double chi2_point = residual.transpose() * point_and_cov.second.inverse() * residual;

        // std::cout << "cov residual \n " <<  point_and_cov.second << std::endl;
        // std::cout << "dist, chi2 " << dist_point <<" , "<<chi2_point << std::endl;

        return std::make_pair(dist_point, chi2_point);
    }

    std::pair<Vector4d, MatrixXd> Track::get_same_invar_pos_and_cov(Vector4d point, double speed_constraint, bool multiple_scattering) const
    {
        (void)speed_constraint;

        double dt = point(this->iv_index) - this->iv_value;
        Vector3d pos_on_track = this->params.segment(0, 3) + this->params.segment(3, 3) * dt;
        Vector4d pos_full = Helper::insertVector(pos_on_track, this->iv_index, point(this->iv_index));

        MatrixXd cov_residual = this->cov.topLeftCorner(3, 3) +
                                this->cov.topRightCorner(3, 3) * dt +
                                this->cov.bottomLeftCorner(3, 3) * dt +
                                this->cov.bottomRightCorner(3, 3) * dt * dt;

        // Multiple scattering matrix. Optional.
        if (multiple_scattering)
        {
            cov_residual += this->Q_block * dt * dt * 4;
        }

        if (this->hit_ids.size() == 4)
        {
            cov_residual *= 4;
        }

        return std::make_pair(pos_full, cov_residual);
    }

    std::pair<double, double> Track::get_same_invar_dist_and_chi2(Vector4d point, double speed_constraint, bool multiple_scattering) const
    {
        (void)speed_constraint;

        auto point_and_cov = get_same_invar_pos_and_cov(point, speed_constraint, multiple_scattering);
        Vector4d residual4 = point_and_cov.first - point;
        auto residual = Helper::removeElement(residual4, this->iv_index);
        double dist_point = residual.norm();
        double chi2_point = residual.transpose() * point_and_cov.second.inverse() * residual;

        // std::cout << "cov \n " <<  point_and_cov.second << std::endl;
        // std::cout << "residual  " <<  residual.transpose() << std::endl;
        // std::cout << "dist, chi2 " << dist_point <<" , "<<chi2_point << std::endl;

        return std::make_pair(dist_point, chi2_point);
    }

    std::pair<Vector4d, MatrixXd> Track::get_cpa_pos_and_cov(Vector4d point, double speed_constraint, bool multiple_scattering) const
    {
        (void)speed_constraint;
        (void)multiple_scattering;

        Vector4d r0 = this->params_full.segment(0, 4);
        Vector4d v0 = this->params_full.segment(4, 4);
        double v0_2 = std::pow(v0.norm(), 2);
        double dt = (point - r0).dot(v0) / v0_2;
        Vector4d closest_point = r0 + v0 * dt;

        MatrixXd J_r0(4, 4), J_v0(4, 4), cov_point(4, 4);

        J_r0.setIdentity();
        J_r0 += -v0 * v0.transpose() / v0_2;

        J_v0.setIdentity();
        J_v0.diagonal() *= dt;
        J_v0 += r0 * (point - r0 - 2 * dt * v0).transpose() / v0_2;

        cov_point = J_r0 * cov_full.topLeftCorner(4, 4) * J_r0.transpose() +
                    J_r0 * cov_full.topRightCorner(4, 4) * J_v0.transpose() +
                    J_v0 * cov_full.bottomLeftCorner(4, 4) * J_r0.transpose() +
                    J_v0 * cov_full.bottomRightCorner(4, 4) * J_v0.transpose();

        return std::make_pair(closest_point, cov_point);
    }

    std::pair<double, double> Track::get_cpa_dist_and_chi2(Vector4d point, double speed_constraint, bool multiple_scattering) const
    {
        auto point_and_cov = get_cpa_pos_and_cov(point, speed_constraint, multiple_scattering);
        Vector4d residual4d = point_and_cov.first - point;
        double dist_point = residual4d.segment(0, 3).norm();
        double chi2_point = residual4d.transpose() * point_and_cov.second.inverse() * residual4d;

        // std::cout << "cov residual \n " << point_and_cov.second << std::endl;
        // std::cout << "dist, chi2: " << dist_point <<" , "<<chi2_point << std::endl;

        return std::make_pair(dist_point, chi2_point);
    }

    // -------------------------------------------------------------------------------------------
    // Class Helper

    // Remove i^th element from a vector and creat a new one.
    Eigen::VectorXd Helper::removeElement(const Eigen::VectorXd &vec4, int i)
    {
        Eigen::VectorXd vec3(vec4.rows() - 1);
        // Ensure the index is valid
        if (i < 0 || i >= vec4.rows())
            throw std::out_of_range("Index i is out of range. Must be between 0 and 3.");

        // Copy elements excluding the i-th element
        for (int j = 0, k = 0; j < vec4.rows(); ++j)
        {
            if (j != i)
                vec3(k++) = vec4(j);
        }
        return vec3;
    }

    // Inseart value to a vector at given index i
    Eigen::VectorXd Helper::insertVector(const Eigen::VectorXd &vec4, int i, float val)
    {
        Eigen::VectorXd vec5(vec4.rows() + 1);

        // Copy elements excluding the i-th element
        for (int j = 0, k = 0; j < vec5.rows(); ++j)
        {
            if (j != i)
                vec5(j) = vec4(k++);
            else
                vec5(j) = val;
        }
        return vec5;
    }

    // Function to insert a row at a specific index
    Eigen::MatrixXd Helper::insertRow(const Eigen::MatrixXd &mat, const Eigen::RowVectorXd &row, int insert_pos)
    {
        // Check if the insert position is valid
        if (insert_pos < 0 || insert_pos > mat.rows())
            throw std::invalid_argument("Invalid insert position, Index i is out of range");

        // Create a new matrix with one additional row
        Eigen::MatrixXd new_mat(mat.rows() + 1, mat.cols());

        new_mat.topRows(insert_pos) = mat.topRows(insert_pos);                                 // Copy the top rows (before the insert position)
        new_mat.row(insert_pos) = row;                                                         // Insert the new row
        new_mat.bottomRows(mat.rows() - insert_pos) = mat.bottomRows(mat.rows() - insert_pos); // Copy the bottom rows (after the insert position)

        return new_mat;
    }

    // Function to insert a column at a specific index
    Eigen::MatrixXd Helper::insertColumn(const Eigen::MatrixXd &mat, const Eigen::VectorXd &col, int insert_pos)
    {
        // Check if the insert position is valid
        if (insert_pos < 0 || insert_pos > mat.cols())
            throw std::invalid_argument("Invalid insert position");

        // Create a new matrix with one additional column
        Eigen::MatrixXd new_mat(mat.rows(), mat.cols() + 1);

        // Copy the left columns (before the insert position)
        new_mat.leftCols(insert_pos) = mat.leftCols(insert_pos);
        // Insert the new column
        new_mat.col(insert_pos) = col;
        // Copy the right columns (after the insert position)
        new_mat.rightCols(mat.cols() - insert_pos) = mat.rightCols(mat.cols() - insert_pos);

        return new_mat;
    }

    // Function to insert a row and column of zeros at the i-th index
    Eigen::MatrixXd Helper::insertRowAndColumnOfZeros(const Eigen::MatrixXd &mat, int i)
    {
        // Check if the insert position is valid
        if (i < 0 || i > mat.rows() || i > mat.cols())
            throw std::invalid_argument("Invalid insert position");

        // Create a new matrix with one additional row and column
        Eigen::MatrixXd new_mat(mat.rows() + 1, mat.cols() + 1);

        // Copy the top-left block (before the i-th row and column)
        new_mat.topLeftCorner(i, i) = mat.topLeftCorner(i, i);
        // Copy the top-right block (before the i-th row, after the i-th column)
        new_mat.topRightCorner(i, mat.cols() - i) = mat.topRightCorner(i, mat.cols() - i);
        // Copy the bottom-left block (after the i-th row, before the i-th column)
        new_mat.bottomLeftCorner(mat.rows() - i, i) = mat.bottomLeftCorner(mat.rows() - i, i);
        // Copy the bottom-right block (after the i-th row and column)
        new_mat.bottomRightCorner(mat.rows() - i, mat.cols() - i) = mat.bottomRightCorner(mat.rows() - i, mat.cols() - i);

        // Set the i-th row and column to zeros
        new_mat.row(i).setZero();
        new_mat.col(i).setZero();

        return new_mat;
    }

    // Class TreeReaderDigi------------------------------------------------------------------------------------------------
    TreeReaderDigi::TreeReaderDigi(std::string filename)
    {
        auto tree_name = "data";
        File = TFile::Open(filename.c_str());
        TreeData = (TTree *)File->Get(tree_name);
        TreeMetadata = (TTree *)File->Get("metadata");

        entries_total = TreeData->GetEntries();

        // Setup tree pointer to data buffer
        TreeData->SetBranchAddress("Digi_x", &Digi_x);
        TreeData->SetBranchAddress("Digi_y", &Digi_y);
        TreeData->SetBranchAddress("Digi_z", &Digi_z);
        TreeData->SetBranchAddress("Digi_t", &Digi_t);
        TreeData->SetBranchAddress("Digi_edep", &Digi_edep);
        TreeData->SetBranchAddress("Digi_trackID", &Digi_trackID);
        TreeData->SetBranchAddress("Digi_pdgID", &Digi_pdgID);
        TreeData->SetBranchAddress("Digi_detectorID", &Digi_detectorID);
        TreeData->SetBranchAddress("Digi_type", &Digi_type);
        TreeData->SetBranchAddress("Digi_hitInds", &Digi_hitInds);
        TreeData->SetBranchAddress("Digi_direction", &Digi_direction);

        // Read metadata
        TreeMetadata->SetBranchAddress("SimulationName", &SimulationName);
        TreeMetadata->SetBranchAddress("Geometry", &Geometry);
        TreeMetadata->SetBranchAddress("Generator", &Generator);
        TreeMetadata->SetBranchAddress("Uncertainty_t", &Uncertainty_t);
        TreeMetadata->SetBranchAddress("Uncertainty_x", &Uncertainty_x);
        TreeMetadata->SetBranchAddress("Uncertainty_y", &Uncertainty_y);
        TreeMetadata->SetBranchAddress("Uncertainty_z", &Uncertainty_z);
        TreeMetadata->GetEntry(0);
        digi_unc << Uncertainty_x, Uncertainty_y, Uncertainty_z, Uncertainty_t;
    }

    std::vector<std::unique_ptr<DigiHit>> TreeReaderDigi::GetEntry(long long entry)
    {
        std::vector<std::unique_ptr<DigiHit>> hits_tmp;

        TreeData->GetEntry(entry);
        for (size_t i = 0; i < (*Digi_x).size(); i++)
        {
            int dir_x = ((*Digi_direction)[i] / 100) % 10;
            int dir_y = ((*Digi_direction)[i] / 10) % 10;
            int dir_z = ((*Digi_direction)[i]) % 10;
            double unc[3] = {0, 0, 0};
            unc[dir_x] = digi_unc[0];
            unc[dir_y] = digi_unc[1];
            unc[dir_z] = digi_unc[2];
            int hit_layer = ((*Digi_detectorID)[i] / 100000) % 1000;
            int hit_detector_group = ((*Digi_detectorID)[i] / 100000000000) % 1000;
            // std::cout<<"Add hit: " << i <<" " << hit_layer << std::endl;
            auto hit = std::make_unique<Tracker::DigiHit>((*Digi_x)[i],
                                                          (*Digi_y)[i],
                                                          (*Digi_z)[i],
                                                          (*Digi_t)[i],
                                                          unc[0], unc[1], unc[2], digi_unc[3],
                                                          dir_z,
                                                          hit_layer,
                                                          hit_detector_group, i);
            hit->detector_id = (*Digi_detectorID)[i];
            hits_tmp.push_back(std::move(hit));
        }
        return hits_tmp;
    }

    std::unordered_map<int, std::vector<DigiHit *>> TreeReaderDigi::ProcessHits(std::vector<std::unique_ptr<DigiHit>> &hits_tmp)
    {
        std::unordered_map<int, std::vector<DigiHit *>> hits_dict;
        // for (auto &hit : hits_tmp)
        // {
        //     auto group = hit->group;
        //     if (group == 0)
        //         continue;
        //     hits_dict[group].push_back(hit.get());
        // }

        std::unordered_map<int, int> hit_inds_to_skip;
        for (int i = 0; i < hits_tmp.size(); i++)
        {
            if (hit_inds_to_skip.count(i) > 0)
            {
                continue;
            }

            // Check if there is an adjacent hit
            bool found_adjacent = false;
            int k = -1;
            for (int j = i + 1; j < hits_tmp.size(); j++)
            {
                if (std::abs(hits_tmp[i]->detector_id - hits_tmp[j]->detector_id) == 1 &&
                    std::abs(hits_tmp[i]->t() - hits_tmp[j]->t()) < 5)
                {
                    found_adjacent = true;
                    k = j;
                    hit_inds_to_skip[k] = 1;
                    break;
                }
            }

            // Alter the position and time of the first hit in the pair
            if (found_adjacent)
            {
                hits_tmp[i]->setx((hits_tmp[i]->x() + hits_tmp[k]->x()) * 0.5);
                hits_tmp[i]->sety((hits_tmp[i]->y() + hits_tmp[k]->y()) * 0.5);
                hits_tmp[i]->setz((hits_tmp[i]->z() + hits_tmp[k]->z()) * 0.5);
                hits_tmp[i]->sett((hits_tmp[i]->t() + hits_tmp[k]->t()) * 0.5);
                hits_tmp[i]->setex(hits_tmp[i]->ex() * 0.707);
                hits_tmp[i]->setey(hits_tmp[i]->ey() * 0.707);
                hits_tmp[i]->setez(hits_tmp[i]->ez() * 0.707);
                hits_tmp[i]->setet(hits_tmp[i]->et() * 0.707);
            }

            auto group = hits_tmp[i]->group;
            if (group == 0)
                continue;
            hits_dict[group].push_back(hits_tmp[i].get());
        }

        return hits_dict;
    }

    // Class TreeWriterRecon------------------------------------------------------------------------------------------------
    TreeWriterRecon::TreeWriterRecon(std::string filename_recon,
                                     std::string filename_digi,
                                     std::string filename_sim,
                                     bool save_raw_reduced) : iroot::file::EntryCopy()
    {

        auto output_tree_name = "data";
        outputFile = TFile::Open(filename_recon.c_str(), "RECREATE");
        outputTree = new TTree(output_tree_name, "Reconstruction Tree");
        outputTreeMetadata = new TTree("metadata", "Metadata for reconstruction");

        // Write metadata
        outputTreeMetadata->Branch("ReconstructionConfigStr", &meta_ReconstructionConfigStr);

        // Setup tree pointer to data buffer
        outputTree->Branch("SimEntry", &SimEntry);
        outputTree->Branch("Track_x0", &Track_x0);
        outputTree->Branch("Track_y0", &Track_y0);
        outputTree->Branch("Track_z0", &Track_z0);
        outputTree->Branch("Track_t0", &Track_t0);
        outputTree->Branch("Track_kx", &Track_kx);
        outputTree->Branch("Track_ky", &Track_ky);
        outputTree->Branch("Track_kz", &Track_kz);
        outputTree->Branch("Track_kt", &Track_kt);
        outputTree->Branch("Track_cov", &Track_cov); // Have to be flattened, each track takes 6x6=36 elements
        outputTree->Branch("Track_chi2", &Track_chi2);
        outputTree->Branch("Track_id", &Track_id);
        outputTree->Branch("Track_iv_ind", &Track_iv_ind);
        outputTree->Branch("Track_iv_err", &Track_iv_err);
        outputTree->Branch("Track_digiInds", &Track_digiInds);
        outputTree->Branch("Vertex_x0", &Vertex_x0);
        outputTree->Branch("Vertex_y0", &Vertex_y0);
        outputTree->Branch("Vertex_z0", &Vertex_z0);
        outputTree->Branch("Vertex_t0", &Vertex_t0);
        outputTree->Branch("Vertex_cov", &Vertex_cov); // Have to be flattened, each track takes 6x6=36 elements
        outputTree->Branch("Vertex_chi2", &Vertex_chi2);
        outputTree->Branch("Vertex_id", &Vertex_id);
        outputTree->Branch("Vertex_trackInds", &Vertex_trackInds);

        // Setup copy methods for raw/digi data
        EN_COPY_DIGI = filename_digi.size() > 0 ? true : false;
        EN_COPY_RAW = filename_sim.size() > 0 ? true : false;
        if (EN_COPY_DIGI)
        {
            digiFile = TFile::Open(filename_digi.c_str());
            digiTree = (TTree *)digiFile->Get(output_tree_name);
            digiTreeMetadata = (TTree *)digiFile->Get("metadata");
            Setup(digiTree, outputTree); // Setup for copying from digiTree to outputTree
        }
        if (EN_COPY_RAW)
        {
            simFile = TFile::Open(filename_sim.c_str());
            simTree = (TTree *)simFile->Get(output_tree_name);
            simTreeMetadata = (TTree *)simFile->Get("metadata");
            SetSimBranches(save_raw_reduced);
            Setup(simTree, outputTree); // Setup for copying from rawTree to outputTree
        }
    }

    void TreeWriterRecon::SetSimBranches(bool save_raw_reduced)
    {
        if (EN_COPY_RAW)
        {
            if (save_raw_reduced)
            {
                std::vector<std::string> branches_enabled = {"Run_number",
                                                             "Evt_number",
                                                             "Evt_weight",
                                                             "Seed_0",
                                                             "Seed_1",
                                                             "Gen_*"};

                // Disable all branches first
                simTree->SetBranchStatus("*", 0);

                // Enable selected ones
                for (auto &br : branches_enabled)
                    simTree->SetBranchStatus(br.c_str(), 1);
            }
        }
    }

    int TreeWriterRecon::ApplyRecon(TrackList &tracks, VertexLilst &vertices, int &simulation_entry_number)
    {
        SimEntry = simulation_entry_number;

        for (auto &track : tracks)
        {
            Track_x0.push_back(track->params_full[0]);
            Track_y0.push_back(track->params_full[1]);
            Track_z0.push_back(track->params_full[2]);
            Track_t0.push_back(track->params_full[3]);
            Track_kx.push_back(track->params_full[4]);
            Track_ky.push_back(track->params_full[5]);
            Track_kz.push_back(track->params_full[6]);
            Track_kt.push_back(track->params_full[7]);
            Track_chi2.push_back(track->chi2);
            Track_id.push_back(track->id);
            Track_iv_ind.push_back(track->iv_index);
            Track_iv_err.push_back(track->iv_error);
            // Flatten cov matrix
            for (Eigen::MatrixXd::Index i = 0; i < track->cov.size(); ++i)
                Track_cov.push_back(track->cov.data()[i]);
            // Flatten hit ids, separate by -1 at the end
            for (auto hitid : track->hit_ids)
                Track_digiInds.push_back(hitid);
            Track_digiInds.push_back(-1);
        }

        for (auto &vertex : vertices)
        {
            Vertex_x0.push_back(vertex->params[0]);
            Vertex_y0.push_back(vertex->params[1]);
            Vertex_z0.push_back(vertex->params[2]);
            Vertex_t0.push_back(vertex->params[3]);
            Vertex_chi2.push_back(vertex->chi2);
            Vertex_id.push_back(vertex->id);
            // Flatten cov matrix
            for (Eigen::MatrixXd::Index i = 0; i < vertex->cov.size(); ++i)
                Vertex_cov.push_back(vertex->cov.data()[i]);
            // Flatten hit ids, separate by -1 at the end
            for (auto trackid : vertex->track_ids)
                Vertex_trackInds.push_back(trackid);
            Vertex_trackInds.push_back(-1);
        }
        return 0;
    }

    int TreeWriterRecon::ApplyCopy(long long entry)
    {
        if (EN_COPY_DIGI)
            ReadSource(digiTree, entry);
        if (EN_COPY_RAW)
            ReadSource(simTree, entry);

        return 0;
    }

    void TreeWriterRecon::Clear()
    {
        Track_x0.clear();
        Track_y0.clear();
        Track_z0.clear();
        Track_t0.clear();
        Track_kx.clear();
        Track_ky.clear();
        Track_kz.clear();
        Track_kt.clear();
        Track_cov.clear();
        Track_chi2.clear();
        Track_id.clear();
        Track_iv_ind.clear();
        Track_iv_err.clear();
        Track_digiInds.clear();
        Vertex_x0.clear();
        Vertex_y0.clear();
        Vertex_z0.clear();
        Vertex_t0.clear();
        Vertex_cov.clear();
        Vertex_chi2.clear();
        Vertex_id.clear();
        Vertex_trackInds.clear();
    }

    void TreeWriterRecon::Fill()
    {
        outputTree->Fill();
        this->Clear();
    }

} // namespace Tracker
