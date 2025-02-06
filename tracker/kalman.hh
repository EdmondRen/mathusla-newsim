#ifndef kalman_HH
#define kalman_HH

#include <iostream>
#include <numeric>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::MatrixXd;
using Eigen::VectorXd;
// typedef Eigen::SparseMatrix<double> MatrixXdSp;
using MatrixXdSp = Eigen::SparseMatrix<double>;

namespace Kalman
{
    class KF_Forward
    {
    public:
        KF_Forward(int ndimMeasure, int ndimStates);

        void SetInitialState(VectorXd &xf0, MatrixXd &Cf0);

        // Update the predicted next state (Xp), state covariance (Cp) and residual covariance (Rp)
        void UpdateMatrix(MatrixXd &Vi, MatrixXdSp &Hi, MatrixXdSp &Fi, MatrixXd &Qi);

        // Calculate the predicted chi2 with the next measurement
        double TryMeasurement(const VectorXd &mi, int n_sigma = -1);

        // Run filter on the next measurement
        double Filter(const VectorXd mi);

        VectorXd GetState() { return xf; }
        MatrixXd GetCov() { return Cf; }
        VectorXd GetPredict() { return m_predict; }
        VectorXd GetPredictErr() { return m_predict_err; }
        double GetChi2Step() { return chi2_step; }
        double GetChi2() { return chi2_total; }

    protected:
        // Dimensions
        int Nmeas, Nstat;
        MatrixXd Imeas, Istat; // residual cov predicted/filtered

        // Measurement
        // Use pointer to avoid memory copy
        VectorXd m_predict;     // predicted Measurement vector
        VectorXd m_predict_err; // predicted Measurement error
        MatrixXdSp *H;          // Measurement matrix H
        MatrixXd *V;            // Measurement uncertainty matrix V

        // State vectors
        VectorXd xp, xf; // xp: predicted, xf: filtered

        // System dynamics
        // Use pointer to avoid memory copy
        MatrixXdSp *F; // Fi: process matrix
        MatrixXd *Q;   // Qi: process uncertainty

        // Uncertainty and residual
        MatrixXd Cp, Cf;         // Variation Matrices. Cp: predicted, Cf: filtered
        VectorXd rp;             // predicted residual
        MatrixXd Rp, Rf, Rp_inv; // residual (of measurement) cov predicted/filtered

        // Chi2
        double chi2_step, chi2_total;
    };

    class KF_FULL
    {
    public:
        // Initialize the filter with
        KF_FULL(int ndimMeasure, int ndimStates);

        void SetInitialState(VectorXd &m0, MatrixXd &V0, VectorXd &xf0, MatrixXd &Cf0, MatrixXd &H0);

        // Run filter on the next measurement
        // Update the predicted next state (Xp), state covariance (Cp) and residual covariance (Rp)
        double FilterStep(VectorXd &mi, MatrixXd &Vi, MatrixXd &Hi, MatrixXd &Fi, MatrixXd &Qi);

        void init_smooth();
        double smooth_step_try();
        void smooth_step(bool drop = false);

        // Run smoothing on all steps without dropping
        void backward_smooth();

        double chift_sum, chism_sum;
        int CURRENT_STEP;
        
        // inline VectorXd GetState() { return xf; }
        // inline MatrixXd GetCov() { return Cf; }
        // inline VectorXd GetPredict() { return m_predict; }
        // inline VectorXd GetPredictErr() { return m_predict_err; }
        // inline double GetChi2Step() { return chi2_step; }
        // inline double GetChi2() { return chi2_total; }

    protected:
        // Dimensions
        int Nmeas, Nstat;
        MatrixXd Imeas, Istat; // residual cov predicted/filtered

        // Suppose there are N+1 measurements starting from index 0
        // measurements
        std::vector<VectorXd> m; // measurements
        std::vector<MatrixXd> V; // measurement uncertainty matrices
        std::vector<MatrixXd> H; // measurement matrices
        // states
        std::vector<VectorXd> xp;  // measurements
        std::vector<VectorXd> xf;  // measurements
        std::vector<VectorXd> xsm; // measurements
        // Extrapolation Matrices
        std::vector<MatrixXd> F; // Extrapolate function              [0...N-1]
        std::vector<MatrixXd> Q; // Extrapolate uncertainty           [0...N-1]
        // Variation Matrices
        std::vector<MatrixXd> Cp;  // state Cov predicted               [_,1...N] *Initial value needs a placeholder
        std::vector<MatrixXd> Cf;  // state Cov filtered                [0...N]
        std::vector<MatrixXd> Csm; // state Cov smoothed                [0...N] *------INIT when smoothing
        // Residuals
        std::vector<VectorXd> rp; // residual predicted                [_,1...N] *Initial value needs a placeholder
        std::vector<MatrixXd> Rp; // residual Cov predicted            [_,1...N] *Initial value needs a placeholder
        // std::vector<MatrixXd> Rf;  // residual Cov filtered             [0...N]
        std::vector<MatrixXd> Rsm; // residual Cov smoothed             [0...N] *------INIT when smoothing
        // Gains
        std::vector<MatrixXd> K; // Forward Kalman gain               [0...N] *Initial value derived from Cf[0] and H[0]
        std::vector<MatrixXd> A; // Backward Kalman gain              [0...N-1]
        // Chi-squares
        std::vector<double> chift; // Forward (filtering) chi2          [0...N] *Initial value is constant, chift[0] = 0
        std::vector<double> chism; // Backward (smoothing) chi2         [0...N] *------INIT when smoothing

        // Temp variable for smoothing
        std::vector<int> STEPS_LIST;

        bool SMOOTH_STEP_TRIED;
        VectorXd xsm_temp;
        MatrixXd Csm_temp;
        MatrixXd Rsm_temp;
        MatrixXd A_temp;
        double chism_temp;
    };

} // namespace kalman

#endif