#ifndef kalman_model_HH
#define kalman_model_HH

#include <iostream>
#include <memory>

#include <Eigen/Dense>
#include <TMinuit.h>

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
        // If the hit passed the range and chi2 cut, return the chi2 value
        // Else, return -1
        float try_measurement(const Tracker::DigiHit &hit, float sigma_cut, float chi2_cut);

        // Add a measurement
        int add_measurement(const Tracker::DigiHit &hit);

        // Update process noise (multiple scattering)
        int update_Q(float step, float multiple_scattering_p = 500, float multiple_scattering_length = 0.06823501107481977);

        // Run filter forward
        std::unique_ptr<Tracker::Track> run_filter(const std::vector<Tracker::DigiHit *> &hits);

        double step_current, step_next, step_size;

    protected:
        // Filter instance
        KF_Forward kf;
        bool DEBUG;

        // Parameters
        int enMultipleScattering;
        int Nmeas, Nstat; // Dimensions

        // Initial state and covariance
        VectorXd xf0; // state, 1x6 col vector
        MatrixXd Cf0; // covariance, 6x6
        // Intermediate parameters
        MatrixXd Vi;   // Measurement uncertaint, 3x3
        MatrixXdSp Hi; // Measurement matrix, 3x6
        MatrixXdSp Fi; // State transfer dynamics, 6x6
        MatrixXd Qi;   // Process noise matrix, 6x6

        // Final output: a track object that holds all the information
        std::unique_ptr<Tracker::Track> track; // state, 1x6 col vector
    };

    class LSVertex4DFitter
    {
    public:
        LSVertex4DFitter() : minimizer(npar),cov(4,4), parameters(npar), parameter_errors(npar)
        {
            minimizer.SetFCN(cost);
            int QUIET_MODE = -1;
            minimizer.SetPrintLevel(QUIET_MODE);
        }
        static const int npar = 4;
        TMinuit minimizer;

        // Cost function to minimize
        // Arguments are required by TMinuit. Only the following two are needed here.
        // - &f: save the function value here
        // - *par: list of parameters
        // It has to be a static member
        static void cost(int &npar, double *gin, double &f, double *par, int iflag);

        // Find vertex that minimize chi2
        bool fit(std::vector<Tracker::Track *> tracks,
                 std::vector<double> arg_guess = {},
                 double tolerance = 0.1,
                 double maxcalls = 50000);

        // Holding fit results
        Vector4d params;
        MatrixXd cov;
        double chi2;
        int status;

        // All tracks to be used. Declared static to be accessible to static cost funtion.
        static std::vector<Tracker::Track *> tracks;

    protected:
        std::vector<double> parameters;
        std::vector<double> parameter_errors;
        double cov_matrix[npar][npar];
    };

    class VertexSeed
    {
    public:
        VertexSeed(Tracker::Track *track1,
                   Tracker::Track *track2) : cov(4, 4), ntracks_found(-1)
        {
            this->tracks = std::make_pair(track1, track2);

            auto guess = Tracker::Track::get_closest_midpoint(*track1, *track2);
            this->distance = guess.first;
            this->vertex_guess = guess.second;
        }

        float GetScore()
        {
            // float c = 299.7; // speed of light [mm/ns]

            // Use LS fit to get the chi2 and covariance matrix of this pair
            auto vertex_fitter = Kalman::LSVertex4DFitter();
            std::vector<Tracker::Track *> tracks_fit = {tracks.first, tracks.second};

            // Reuse the calculated best guess to initialize the fitter
            std::vector<double> vertex_guess_vec(4);
            std::copy(vertex_guess.data(), vertex_guess.data() + vertex_guess.size(), vertex_guess_vec.begin());
            vertex_fitter.fit(tracks_fit, vertex_guess_vec, 0.1, 1000);

            // Get vertex
            this->vertex_fit = vertex_fitter.params;

            // Get chi2
            this->chi2 = vertex_fitter.chi2;

            // Get covariance
            this->cov = vertex_fitter.cov;

            // Make a score out of it
            this->score = -chi2;

            return this->score;
        }

        // Required data
        std::pair<Tracker::Track *, Tracker::Track *> tracks;
        MatrixXd cov;
        Eigen::Vector4d vertex_guess, vertex_fit;
        float score;
        float distance, chi2; // Distance between the two tracks
        int ntracks_found;    // Number of track compatible with this seed
    };

    class KalmanVertex4D
    {
    public:
        // Constructor
        KalmanVertex4D(bool debug = false);

        // Initialize
        int init_state(VertexSeed &seed);

        // Try a measurement
        // If the hit passed the range and chi2 cut, return the chi2 value
        // Else, return -1
        float try_measurement(const Tracker::Track &track, float distance_cut, float chi2_cut);

        // Add a measurement
        int add_measurement();

        // Update process noise (multiple scattering)
        int update_Q(float step, float multiple_scattering_p = 500, float multiple_scattering_length = 0.06823501107481977);

        // Run filter forward
        std::unique_ptr<Tracker::Vertex> run_filter(const std::vector<Tracker::Track *> &tracks, VertexSeed *seed = nullptr);

        double step_current, step_next, step_size;

    protected:
        // Filter instance
        KF_Forward kf;
        bool DEBUG;

        // Parameters
        int enMultipleScattering;
        int Nmeas, Nstat; // Dimensions

        // Initial state and covariance
        VectorXd mi;  // measurement, 1x4 col vector
        VectorXd xf0; // state, 1x4 col vector
        MatrixXd Cf0; // covariance, 4x4
        // Intermediate parameters
        MatrixXd Vi;   // Measurement uncertaint, 4x4
        MatrixXdSp Hi; // Measurement matrix, 4x4
        MatrixXdSp Fi; // State transfer dynamics, 4x4
        MatrixXd Qi;   // Process noise matrix, 4x4

        // Final output: a vertex object that holds all the information
        std::unique_ptr<Tracker::Vertex> vertex; // state, 1x6 col vector
    };

} // namespace Kalman

#endif
