#include "types.hh"

namespace Tracker
{
    // -------------------------------------------------------------------------------------------
    // Class Track

    VectorXd Track::convert_to_time()
    {
        if (param_ind < 3)
        {
            params_time = VectorXd(6);
            cov_time = MatrixXd(6, 6);
            t0 = params(2);
            MatrixXdSp jac(6, 6);

            // Setting parameters
            for (int j = 0, k = 0; j < 3; ++j)
            {
                if (j == param_ind)
                {
                    params_time(j) = param_value;
                    params_time(j + 3) = 1 / params(5);
                }
                else
                {
                    params_time(j) = params(k++);
                    params_time(j + 3) = params(j + 3) / params(5);
                }
            }

            // Setting Jacobian
            for (int j = 0; j < 3; ++j)
            {
                if (j == param_ind)
                {
                    jac.insert(j, j) = params_time(j + 3);
                    jac.insert(j + 3, j + 3) = -params_time(j + 3) * params_time(j + 3);
                }
                else
                {
                    jac.insert(j, j) = 1;
                    jac.insert(j, param_ind) = params_time(j + 3);
                    jac.insert(j + 3, j + 3) = params_time(param_ind + 3);
                    jac.insert(j + 3, param_ind + 3) = -params_time(j + 3) * params_time(param_ind + 3);
                }
            }
            cov_time = jac * cov * jac.transpose();
        }
        else
        {
            params_time = params; // Do nothing if it is already using time as independent variable.
            cov_time = cov;
            t0 = param_value;
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

    std::pair<double, double> Track::get_closest_point_dist_chi2(Vector4d point, double speed_constraint)
    {
        double dt = point(3) - this->t0;
        Vector3d pos_new = this->params_time.segment(0, 3) + this->params_time.segment(3, 3) * dt;
        Vector3d pos_residual = pos_new - point.segment(0, 3);
        MatrixXd cov_residual = this->cov_time.topLeftCorner(3, 3) +
                                this->cov_time.topRightCorner(3, 3) * dt +
                                this->cov_time.bottomLeftCorner(3, 3) * dt +
                                this->cov_time.bottomRightCorner(3, 3) * dt * dt;
        double chi2_point = pos_residual.transpose() * cov_residual.inverse() * pos_residual;
        double dist_point = pos_residual.norm();

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

} // namespace Tracker
