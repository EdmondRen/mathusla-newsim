#ifndef tracker_types_HH
#define tracker_types_HH

#include <iostream>
#include <memory>

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
        static Eigen::VectorXd removeElement(const Eigen::VectorXd &vec4, int i);
        // Inseart value to a vector at given index i
        static Eigen::VectorXd insertVector(const Eigen::VectorXd &vec4, int i, float val);
        // Function to insert a row at a specific index
        static Eigen::MatrixXd insertRow(const Eigen::MatrixXd &mat, const Eigen::RowVectorXd &row, int insert_pos);
        // Function to insert a column at a specific index
        static Eigen::MatrixXd insertColumn(const Eigen::MatrixXd &mat, const Eigen::VectorXd &col, int insert_pos);
        // Function to insert a row and column of zeros at the i-th index
        static Eigen::MatrixXd insertRowAndColumnOfZeros(const Eigen::MatrixXd &mat, int i);
    };

    class DigiHit
    {
    public:
        DigiHit(double _x, double _y, double _z, double _t, double _ex, double _ey, double _ez, double _et, int _iv_index, int _layer, int _group) : iv_index(_iv_index), layer(_layer), group(_group)
        {
            vec4 << _x, _y, _z, _t;
            vec4_err << _ex, _ey, _ez, _et;
            update_vec3();
        }

        // Update 3-vector
        // Must be called after changing the vec4 value
        void update_vec3()
        {
            vec3 = Helper::removeElement(vec4, iv_index);
            vec3_err = Helper::removeElement(vec4_err, iv_index);
        }

        // Required data
        Vector4d vec4, vec4_err; // {x,y,z,t} and their uncertainty
        Vector3d vec3, vec3_err; // Independent variable removed
        int iv_index;           // Index of the independent variable. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int layer;               // Which detector layer this hit belongs to
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
        inline double get_step() const { return vec4(iv_index); }
        inline double get_steperr() const { return vec4_err(iv_index); }
        inline Vector3d get_vec3() const { return vec3; }
        inline Vector3d get_err3() const { return vec3_err; }
    };
    using HitList = std::vector<std::unique_ptr<DigiHit>>;

    class Track
    {
    public:
        Track() {}

        // Required data
        VectorXd params; // {x0, y0, z0, t0, Ax, Ay, Az, At}, except for the ones that is used as independent variable
        MatrixXd cov;    // covariance of {x0, y0, z0, t0, Ax, Ay, Az, At}, except for the ones that is used as independent variable
        float chi2;
        int iv_index;   // Index of the independent variable. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int iv_value; // Value of the independent variable
        int iv_error; // Uncertainty of the independent variable

        // Optional data
        int id;
        std::vector<int> hit_ids;

        // Parameters based on time
        float t0, t0_err;
        VectorXd params_time; // {x0, y0, z0, vx, vy, vz}
        MatrixXd cov_time;    // covariance of {x0, y0, z0, vx, vy, vz}

        // Full 8-D parameters
        VectorXd params_full; // {x0, y0, z0, t0, Ax, Ay, Az, At}
        MatrixXd cov_full;    // covariance of {x0, y0, z0, t0, Ax, Ay, Az, At}        

        // Convert to using time as independent variable
        VectorXd convert_to_time();

        // Calculate closest-approach midpoint and distance
        // Return: pair of <closest distance, {x,y,z,t}>
        static std::pair<double, Vector4d> get_closest_midpoint(const Track &track1, const Track &track2);

        // Get the <dist, chi2> at the same time of a given point
        std::pair<double, double> get_same_time_dist_and_chi2(Vector4d point, double speed_constraint = -1.0) const;
    
        // Get the <dist, cov> at the closest point of approach (CPA) on the track to a given point
        std::pair<Vector4d, MatrixXd> get_cpa_pos_and_cov(Vector4d point, double speed_constraint = -1.0) const;        

        // Get the <dist, chi2> at the closest point of approach (CPA) on the track to a given point
        std::pair<double, double> get_cpa_dist_and_chi2(Vector4d point, double speed_constraint = -1.0) const;        
    };
    using TrackList = std::vector<std::unique_ptr<Track>>;

    class TrackSeed
    {
    public:
        TrackSeed(DigiHit *hit1, DigiHit *hit2) : nhits_found(-1)
        {
            float c = 299.7; // speed of light [mm/ns]
            VectorXd dvec = hit2->vec4.array() - hit1->vec4.array();
            this->dstep = std::abs(hit2->get_step() - hit1->get_step());
            this->dr = dvec.segment(0, 3).norm();
            this->dt = std::abs(dvec(3));
            this->score = std::abs(dr / c - dt);

            // Make a pair and put the earlier one in front
            this->hits = hit1->t() < hit2->t() ? std::make_pair(hit1, hit2) : std::make_pair(hit2, hit1);
        }

        // Required data
        std::pair<DigiHit *, DigiHit *> hits;
        float score;
        float dr, dt, dstep;
        int nhits_found; // Number of hits found using this seed. Save time when reusing this seed.

        // float GetScore() { return score; }
    };

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
        std::vector<Track *> tracks;
    };
    using VertexLilst = std::vector<std::unique_ptr<Vertex>>;

} // namespace Tracker

#endif
