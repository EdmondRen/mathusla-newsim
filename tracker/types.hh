#ifndef tracker_types_HH
#define tracker_types_HH

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::VectorXd;
using MatrixXdSp = Eigen::SparseMatrix<double>;

namespace Tracker
{

    // Helper functions
    class Helper
    {
    public:
        Helper() {}
        // Remove i^th element from a vector and creat a new one.
        static Eigen::VectorXd removeElement(const Eigen::VectorXd &vec4, int i)
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
        static Eigen::VectorXd insertVector(const Eigen::VectorXd &vec4, int i, float val)
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
        static Eigen::MatrixXd insertRow(const Eigen::MatrixXd &mat, const Eigen::RowVectorXd &row, int insert_pos)
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
        static Eigen::MatrixXd insertColumn(const Eigen::MatrixXd &mat, const Eigen::VectorXd &col, int insert_pos)
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
        static Eigen::MatrixXd insertRowAndColumnOfZeros(const Eigen::MatrixXd &mat, int i)
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
    };

    class DigiHit
    {
    public:
        DigiHit(double _x, double _y, double _z, double _t, double _ex, double _ey, double _ez, double _et, int _param_ind, int _group) : param_ind(_param_ind), group(_group)
        {
            vec4 << _x, _y, _z, _t;
            vec4_err << _ex, _ey, _ez, _et;
            update_vec3();
        }

        // Update 3-vector
        // Must be called after changing the vec4 value
        void update_vec3()
        {
            vec3 = Helper::removeElement(vec4, param_ind);
            vec3_err = Helper::removeElement(vec4_err, param_ind);
        }

        // Required data
        Vector4d vec4, vec4_err; // {x,y,z,t} and their uncertainty
        Vector3d vec3, vec3_err; // Independent variable removed
        int param_ind;           // Index of the independent variable. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int group;               // Which detector grop this hit belongs to

        // Optional
        int id;
        int type;
        u_int64_t detector_id;

        // Set and Get methods
        inline double x() const { return vec4(0); }
        inline double y() const { return vec4(1); }
        inline double z() const { return vec4(2); }
        inline double t() const { return vec4(3); }
        inline void setx(float x) { vec4(0) = x; }
        inline void sety(float x) { vec4(1) = x; }
        inline void setz(float x) { vec4(2) = x; }
        inline void sett(float x) { vec4(3) = x; }
        inline double ex() const { return vec4_err(0); }
        inline double ey() const { return vec4_err(1); }
        inline double ez() const { return vec4_err(2); }
        inline double et() const { return vec4_err(3); }
        inline void setex(float x) { vec4_err(0) = x; }
        inline void setey(float x) { vec4_err(1) = x; }
        inline void setez(float x) { vec4_err(2) = x; }
        inline void setet(float x) { vec4_err(3) = x; }
        inline double get_step() const { return vec4(param_ind); }
        inline double get_steperr() const { return vec4_err(param_ind); }
        inline Vector3d get_vec3() const { return vec3; }
        inline Vector3d get_err3() const { return vec3_err; }
    };
    using HitList = std::vector<std::unique_ptr<DigiHit>>;

    class Track
    {
    public:
        Track() {}

        // Required data
        VectorXd params; // {x0, y0, z0, t0, Ax, Ay, Az, At}, except for the one that is used as independent parameter
        MatrixXd cov;    // covariance of {x0, y0, z0, t0, Ax, Ay, Az, At}, except for the one that is used as independent parameter
        float chi2;
        int param_ind; // Index of the independent parameter. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.

        // Optional data
        int id;
        std::vector<int> hit_ids;

        // Methods

        // Get the closest point of approach (CPA) of the track to a given point
        double get_closest_point(Vector4d point)
        {

        }
    };
    using TrackList = std::vector<std::unique_ptr<Track>>;

    class Vertex
    {
    public:
        Vertex() {}

        // Required data
        VectorXd params; // {x0, y0, z0, t0}
        MatrixXd cov;    // covariance of {x0, y0, z0, t0}
        float chi2;

        // Optional
        int index;
        std::vector<std::unique_ptr<Track>> tracks;
    };

} // namespace Tracker

#endif
