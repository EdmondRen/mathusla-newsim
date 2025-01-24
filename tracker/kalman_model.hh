#ifndef kalman_model_HH
#define kalman_model_HH

#include <iostream>

#include "Eigen/Dense"

#include "types.hh"
#include "kalman.hh"

namespace Kalman
{

    class KalmanTrack4D
    {
    public:
        // Initialize
        //  multiple_scattering: 0 = OFF, 1 = ON
        KalmanTrack4D(int multiple_scattering,
                      int ndimMeasure = 4,
                      int ndimStates = 6,
                      bool debug = false);

        // Initialize
        int init_state(const Tracker::DigiHit &hit1, const Tracker::DigiHit &hit2, bool use_first = true);

        // make a new step. This will only update the related matrix without further calucalation.
        int new_step(const Tracker::DigiHit &hit);

        // Try a measurement
        int try_measurement(const Tracker::DigiHit &hit, float sigma_cut, float chi2_cut);

        // Add a measurement
        int add_measurement(const Tracker::DigiHit &hit);

        // Update process noise (multiple scattering)
        int update_Q(float step, float multiple_scattering_p = 500, float multiple_scattering_length = 0.06823501107481977);

        // Run filter forward
        Tracker::Track *run_filter(Tracker::HitList hits);

    protected:
        // Filter instance
        KF_Forward kf;

        // Parameters
        int enMultipleScattering;
        int Nmeas, Nstat; // Dimensions

        // Initial state and covariance
        VectorXd xf0; // state, 1x6 col vector
        MatrixXd Cf0; // covariance, 6x6
        // Intermediate parameters
        double step_current, step_next, step_size;
        MatrixXd Vi; // Measurement uncertaint, 3x3
        MatrixXdSp Hi; // Measurement matrix, 3x6
        MatrixXdSp Fi; // State transfer dynamics, 6x6
        MatrixXd Qi; // Process noise matrix, 6x6

        // Final output: a track object that holds all the information
        Tracker::Track *track_recon; // state, 1x6 col vector

        bool DEBUG;
    };

} // namespace Kalman

#endif