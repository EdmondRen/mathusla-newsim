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
        // Initialize the filter with
        KF_Forward(int ndimMeasure,
                   int ndimStates) : Nmeas(ndimMeasure),
                                     Nstat(ndimStates),
                                     Imeas(Nmeas, Nmeas),
                                     Istat(Nstat, Nstat),
                                     xp(Nstat),
                                     Cp(Nstat, Nstat),
                                     rp(Nstat),
                                     Rp(Nmeas, Nmeas),
                                     Rf(Nmeas, Nmeas),
                                     Rp_inv(Nmeas, Nmeas),
                                     chi2_step(0)
        {
            // Make some identity matrix
            Imeas.setIdentity();
            Istat.setIdentity();
        }

        void SetInitialState(VectorXd &xf0, MatrixXd &Cf0)
        {
            xf = xf0;
            Cf = Cf0;
            chi2_total = 0;
            chi2_step = 0;
        }

        // Update the predicted next state (Xp), state covariance (Cp) and residual covariance (Rp)
        void UpdateMatrix(MatrixXd &Vi, MatrixXdSp &Hi, MatrixXdSp &Fi, MatrixXd &Qi)
        {
            this->V = &Vi;
            this->H = &Hi;
            this->Q = &Qi;

            this->xp = Fi * this->xf;
            this->Cp = Fi * this->Cf * Fi.transpose() + Qi;
            this->Rp = Vi + Hi * Cp * Hi.transpose();
            this->Rp_inv = this->Rp.inverse();
            this->m_predict = Hi * this->xp;
            this->m_predict_err = Rp.diagonal().array().sqrt();
        }

        // Calculate the predicted chi2 with the next measurement
        double TryMeasurement(const VectorXd &mi, int n_sigma = -1)
        {
            auto rp_i = mi - this->m_predict;

            // Return a very large value if it is obviously out of range
            // This is simply for speed consideration, it
            // avoids preceeding to the chi2 calculation.
            if (n_sigma > 0)
            {
                if (std::abs(rp_i(0)) > n_sigma * m_predict_err(0) ||
                    std::abs(rp_i(1)) > n_sigma * m_predict_err(1) ||
                    std::abs(rp_i(2)) > n_sigma * m_predict_err(2))
                    return 999999;
            }

            double chi2_predict = rp_i.transpose() * this->Rp_inv * rp_i;
            return chi2_predict;
        }

        // Run filter on the next measurement
        double Filter(const VectorXd mi)
        {
            // Get residual
            auto rp_new = mi - (*this->H) * this->xp;

            // Kalman Gain matrix K
            auto K = this->Cp * (*this->H).transpose() * this->Rp_inv;

            // Filtered State
            this->xf = this->xp + K * rp_new; // Combination of the predicted state, measured values, covariance matrix and Kalman Gain

            // Filtered Covariance Matrix
            this->Cf = (this->Istat - K * (*this->H)) * this->Cp;

            // Filtered residual (rf) and residual Covariance matrix (Rf)
            auto rf = mi - (*this->H) * this->xf;
            this->Rf = (this->Imeas - (*this->H) * K) * (*this->V);

            // Chi2 contribution
            this->chi2_step = rf.transpose() * this->Rf.inverse() * rf;
            this->chi2_total += chi2_step;

            return this->chi2_step;
        }

        inline VectorXd GetState() { return xf; }
        inline MatrixXd GetCov() { return Cf; }
        inline VectorXd GetPredict() { return m_predict; }
        inline VectorXd GetPredictErr() { return m_predict_err; }
        inline double GetChi2Step() { return chi2_step; }
        inline double GetChi2() { return chi2_total; }

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
        KF_FULL(int ndimMeasure,
                int ndimStates) : Nmeas(ndimMeasure),
                                  Nstat(ndimStates),
                                  Imeas(Nmeas, Nmeas),
                                  Istat(Nstat, Nstat)
        {
            // Make some identity matrix
            Imeas.setIdentity();
            Istat.setIdentity();

            chift_sum = 0;
            chism_sum = 0;
        }

        void SetInitialState(VectorXd &m0, MatrixXd &V0, VectorXd &xf0, MatrixXd &Cf0, MatrixXd &H0)
        {
            this->m.push_back(m0);
            this->V.push_back(V0);
            this->H.push_back(H0);
            this->xf.push_back(xf0);
            this->Cf.push_back(Cf0);
            // this->Rf.push_back(Rf0);

            auto K0 = Cf0 * H0.transpose() * V0.inverse();
            this->K.push_back(K0);
            this->chift.push_back(0);

            // For those "predicted" variables, fill them with initial state to align the index
            this->xp.push_back(xf0);
            this->Cp.push_back(Cf0);
            this->rp.push_back(m0 * 0.);
        }

        // Run filter on the next measurement
        // Update the predicted next state (Xp), state covariance (Cp) and residual covariance (Rp)
        double FilterStep(VectorXd &mi, MatrixXd &Vi, MatrixXd &Hi, MatrixXd &Fi, MatrixXd &Qi)
        {
            // Append input parameters
            this->m.push_back(mi);
            this->V.push_back(Vi);
            this->H.push_back(Hi);
            this->F.push_back(Fi);
            this->Q.push_back(Qi);

            auto xp_i = Fi * this->xf.back();
            auto Cp_i = Fi * this->Cf.back() * Fi.transpose() + Qi;
            auto Rp_i = Vi + Hi * Cp_i * Hi.transpose();
            auto Rp_inv_i = Rp_i.inverse();
            auto m_predict_i = Hi * xp_i;
            auto m_predict_err_i = Rp_i.diagonal().array().sqrt();

            // Residual
            auto rp_i = mi - m_predict_i;
            double chi2_predict = rp_i.transpose() * Rp_inv_i * rp_i;

            // Kalman Gain matrix K
            auto K_i = Cp_i * Hi.transpose() * Rp_inv_i;

            // Filtered State
            auto xf_i = xp_i + K_i * rp_i; // Combination of the predicted state, measured values, covariance matrix and Kalman Gain

            // Filtered Covariance Matrix
            auto Cf_i = (this->Istat - K_i * Hi) * Cp_i;

            // Filtered residual (rf) and residual Covariance matrix (Rf)
            auto rf_i = mi - Hi * xf_i;
            auto Rf_i = (this->Imeas - Hi * K_i) * Vi;

            // Chi2 contribution
            double chi2_step = rf_i.transpose() * Rf_i.inverse() * rf_i;
            // Save matrices
            this->xp.push_back(xp_i);
            this->Cp.push_back(Cp_i);
            this->K.push_back(K_i);
            this->xf.push_back(xf_i);
            this->Cf.push_back(Cf_i);
            // this->Rf.push_back(Rf_i);
            this->chift.push_back(chi2_step);
            chift_sum += chi2_step;

            return chi2_step;
        }

        // Initialize smoothing
        void init_smooth()
        {
            // Initialized with the last filtered state (X) and Covariance (C)
            xsm.push_back(xf.back());
            Csm.push_back(Cf.back());
            Rsm.push_back(MatrixXd::Zero(Nmeas, Nmeas));
            chism.push_back(0);

            // Make a list from N-1 to 0 (inclusive)
            this->STEPS_LIST = std::vector<int>(xf.size() - 2);
            std::iota(STEPS_LIST.rbegin(), STEPS_LIST.rend(), 0);
            CURRENT_STEP = STEPS_LIST[0];
        }

        double smooth_step_try()
        {
            int i = CURRENT_STEP;
            if (i == -1)
            {
                std::cout << "  KF::Smoothing done, the current state is already the first step" << std::endl;
                return -1;
            }

            // Kalman Gain A
            try
            {
                A_temp = Cf[i] * F[i].transpose() * Cp[i + 1].inverse();
            }
            catch (std::exception &e)
            {
                // Use diagonal approximation if inversion fails
                A_temp = Cf[i] * F[i].transpose() * Cp[i + 1].diagonal().asDiagonal().inverse();
                std::cout << "  KF: Error during smoothing: " << e.what() << std::endl;
            }

            // State
            xsm_temp = xf[i] + A_temp * (this->xsm[0] - xp[i + 1]);
            // Covariance Matrix
            Csm_temp = Cf[i] + A_temp * (this->Csm[0] - Cp[i + 1]) * A_temp.transpose();
            // Residual and Cov of Residual
            VectorXd rsm_temp = m[i] - H[i] * xsm_temp;
            Rsm_temp = V[i] - H[i] * Csm_temp * H[i].transpose();
            // Chi2 contribution
            chism_temp = rsm_temp.transpose() * Rsm_temp.inverse() * rsm_temp;

            // Record the result to be used again in "smooth_step()"
            this->SMOOTH_STEP_TRIED = true;

            return chism_temp;
        }

        void smooth_step(bool drop = false)
        {
            int i = CURRENT_STEP;
            if (i == -1)
            {
                std::cout << "  KF: Smoothing done, the current state is already the first step" << std::endl;
                return;
            }

            // Use the other function to run smoothing step
            if (!SMOOTH_STEP_TRIED)
            {
                smooth_step_try();
            }
            SMOOTH_STEP_TRIED = false;

            // If user choose to drop this step,
            //  revert the contribution to the smoothed state
            if (drop)
            {
                // Modified Kalman Gain K (equation (12b) in Fruhwirth )
                MatrixXd Kmod = Csm_temp * H[i].transpose() * (-V[i] + H[i] * Csm_temp * H[i].transpose()).inverse();
                // Change the smoothed state
                xsm_temp = xsm_temp + Kmod * (m[i] - H[i] * xsm_temp);
                // Change the smoothed Covariance
                Csm_temp = (Istat - Kmod * H[i]) * Csm_temp;
                chism_temp = 0;
            }

            // Insert smoothed values
            xsm.insert(xsm.begin(), xsm_temp);
            Csm.insert(Csm.begin(), Csm_temp);
            Rsm.insert(Rsm.begin(), Rsm_temp);
            A.insert(A.begin(), A_temp);
            chism.insert(chism.begin(), chism_temp);

            // Total chi-square
            chism_sum = accumulate(chism.begin(), chism.end(), 0.0);

            // Move backward by one step
            CURRENT_STEP--;
        }

        void backward_smooth()
        {
            init_smooth();
            while (CURRENT_STEP >= 0)
            {
                smooth_step();
            }
        }

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
        std::vector<VectorXd> rp;  // residual predicted                [_,1...N] *Initial value needs a placeholder
        std::vector<MatrixXd> Rp;  // residual Cov predicted            [_,1...N] *Initial value needs a placeholder
        // std::vector<MatrixXd> Rf;  // residual Cov filtered             [0...N]
        std::vector<MatrixXd> Rsm; // residual Cov smoothed             [0...N] *------INIT when smoothing
        // Gains
        std::vector<MatrixXd> K; // Forward Kalman gain               [0...N] *Initial value derived from Cf[0] and H[0]
        std::vector<MatrixXd> A; // Backward Kalman gain              [0...N-1]
        // Chi-squares
        std::vector<double> chift; // Forward (filtering) chi2          [0...N] *Initial value is constant, chift[0] = 0
        std::vector<double> chism; // Backward (smoothing) chi2         [0...N] *------INIT when smoothing

    public:
        double chift_sum, chism_sum;

        // Temp variable for smoothing
        std::vector<int> STEPS_LIST;
        int CURRENT_STEP;
        bool SMOOTH_STEP_TRIED;

        VectorXd xsm_temp;
        MatrixXd Csm_temp;
        MatrixXd Rsm_temp;
        MatrixXd A_temp;
        double chism_temp;
    };

} // namespace kalman

#endif