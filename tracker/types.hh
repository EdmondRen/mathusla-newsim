#ifndef tracker_types_HH
#define tracker_types_HH

#include <iostream>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::VectorXd;

namespace Matrix
{   

} // namespace Matrix



namespace Tracker
{

    class DigiHit
    {
    public:
        DigiHit(double _x, double _y, double _z, double _t, double _ex, double _ey, double _ez, double _et, int _param_ind, int _group) : param_ind(_param_ind), group(_group)
        {
            vec4 << _x, _y, _z, _t;
            vec4_err << _ex, _ey, _ez, _et;
            update_vec3();
        }

        // Remove an element from a vector and creat a new one.
        Eigen::VectorXd removeElement(const Eigen::VectorXd &vec4, int i)
        {
            Eigen::VectorXd vec3(vec4.rows()-1);

            // Ensure the index is valid
            if (i < 0 || i >= vec4.rows())
            {
                throw std::out_of_range("Index i is out of range. Must be between 0 and 3.");
            }

            // Copy elements excluding the i-th element
            for (int j = 0, k = 0; j < vec4.rows(); ++j)
            {
                if (j != i)
                {
                    vec3(k++) = vec4(j);
                }
            }

            return vec3;
        }        

        // Update 3-vector
        // Must be called after changing the vec4 value
        void update_vec3()
        {
            vec3 = removeElement(vec4, param_ind);
            vec3_err = removeElement(vec4_err, param_ind);
        }

        // Required data
        Vector4d vec4, vec4_err; // {x,y,z,t} and their uncertainty
        Vector3d vec3, vec3_err; // Independent variable removed
        int param_ind;           // Index of the independent variable. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int group; // Which detector grop this hit belongs to

        // Optional
        int index;
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
        inline Vector3d get_vec3() const { return vec3;}
        inline Vector3d get_err3() const { return vec3_err;}
    };
    using HitList = std::vector<DigiHit *>;

    class Track
    {
    public:
        Track() {}

        // Required data
        VectorXd params; // {x0, y0, z0, t0, Ax, Ay, Az, At}
        MatrixXd cov;    // covariance of {x0, y0, z0, t0, Ax, Ay, Az, At}
        float chi2;
        int param_ind; // Index of the independent parameter. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.

        // Optional
        int index;
        std::vector<DigiHit *> hits;
    };
    using TrackList = std::vector<Track *>;

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

    class TrackSeed
    {
    public:
        TrackSeed() {}

        // Required data
        std::pair<DigiHit *, DigiHit *> hits;
        float score;

        float Score() {return 0;}
        float GetScore() { return score; }
    };

} // namespace Tracker

#endif
