#ifndef kalman_HH
#define kalman_HH

#include <iostream>

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
        double TryMeasurement(const VectorXd &mi, int n_sigma = 5)
        {
            auto rp_i = mi - this->m_predict;

            // Return a very large value if it is obviously out of range
            // This is simply for speed consideration, it
            // avoids preceeding to the chi2 calculation.
            if (std::abs(rp_i(0)) > n_sigma * m_predict_err(0) ||
                std::abs(rp_i(1)) > n_sigma * m_predict_err(1) ||
                std::abs(rp_i(2)) > n_sigma * m_predict_err(2))
                return 999999;

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

} // namespace kalman

#endif