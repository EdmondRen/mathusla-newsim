#ifndef tracker_types_HH
#define tracker_types_HH

#include <iostream>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::Vector4d;
using Eigen::VectorXd;

namespace Tracker
{

    class DigiHit
    {
    public:
        DigiHit() {}
        DigiHit(double _x, double _y, double _z, double _t, double _ex, double _ey, double _ez, double _et, int _param_ind, int _group) : param_ind(_param_ind), group(_group)
        {
            vec4 << _x, _y, _z, _t;
            vec4_err << _ex, _ey, _ez, _et;
        }

        // Required data
        Vector4d vec4, vec4_err; // Pack {x,y,z,t} into a vector
        int param_ind;           // Index of the independent parameter. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int group;

        // Optional
        int index;
        int type;
        u_int64_t detector_id;

        // Set and Get methods
        inline double x() { return vec4(0); }
        inline double y() { return vec4(1); }
        inline double z() { return vec4(2); }
        inline double t() { return vec4(3); }
        inline void setx(float x) { vec4(0) = x; }
        inline void sety(float x) { vec4(1) = x; }
        inline void setz(float x) { vec4(2) = x; }
        inline void sett(float x) { vec4(3) = x; }
        inline double ex() { return vec4_err(0); }
        inline double ey() { return vec4_err(1); }
        inline double ez() { return vec4_err(2); }
        inline double et() { return vec4_err(3); }
        inline void setex(float x) { vec4_err(0) = x; }
        inline void setey(float x) { vec4_err(1) = x; }
        inline void setez(float x) { vec4_err(2) = x; }
        inline void setet(float x) { vec4_err(3) = x; }
        inline double get_step() { return vec4(param_ind); }
    };

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

        float Score() {}
        float GetScore() { return score; }
    };

} // namespace Tracker

#endif
