#ifndef kalman_HH
#define kalman_HH

#include <iostream>

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace Kalman
{
    class KF_Forward
    {
    public:
        // Initialize the filter with
        KF_Forward(VectorXd &xf0,
                   MatrixXd &Cf0,
                   int ndimMeasure) : Nmeas(ndimMeasure),
                                      Nstat(xf0.rows()),
                                      Imeas(Nmeas, Nmeas),
                                      Istat(Nstat, Nstat),
                                      xf(xf0),
                                      Cf(Cf0),
                                      xp(Nstat),
                                      Cp(Nstat, Nstat),
                                      rp(Nstat),
                                      Rp(Nstat, Nstat),
                                      Rf(Nstat, Nstat),
                                      Rp_inv(Nstat, Nstat),
                                      chi2_step(0)
        {
            // Make some identity matrix
            Imeas.setIdentity();
            Istat.setIdentity();
        }

        // Update the predicted next state (Xp), state covariance (Cp) and residual covariance (Rp)
        void UpdateMatrix(MatrixXd &Vi, MatrixXd &Hi, MatrixXd &Fi, MatrixXd &Qi)
        {
            this->xp = Fi * this->xf;
            this->Cp = Fi * this->Cf * Fi.transpose() + Qi;
            this->Rp = Vi + Hi * Cp * Hi.transpose();
            this->Rp_inv = this->Rp.inverse();

            this->V = &Vi;
            this->H = &Hi;
            this->Q = &Qi;
        }

        // Calculate the predicted chi2 with the next measurement
        double GetChi2(VectorXd &mi)
        {
            auto rp_i = mi - *this->H * this->xp;
            double chi2_predict = rp_i.transpose() * this->Rp_inv * rp_i;
            return chi2_predict;
        }

        // Run filter on the next measurement
        double Filter(VectorXd &mi)
        {
            // Get residual
            auto rp_new = mi - (*this->H) * this->xp;

            // Kalman Gain matrix K
            auto K = this->Cp * (*this->H).transpose() * this->Rp_inv;

            // Filtered State
            auto xf_new = this->xp + K * rp_new; // Combination of the predicted state, measured values, covariance matrix and Kalman Gain

            // Filtered Covariance Matrix
            this->Cf = (this->Istat - K * (*this->H)) * this->Cp;

            // Filtered residual (rf) and residual Covariance matrix
            auto rf = mi - (*this->H) * xf;
            this->Rf = (this->Imeas - (*this->H) * K) * (*this->V);

            // Chi2 contribution
            this->chi2_step = rf.transpose() * this->Rf.inverse() * rf;
            this->chi2_total += chi2_step;
        }

    protected:
        // Dimensions
        int Nmeas, Nstat;
        MatrixXd Imeas, Istat; // residual cov predicted/filtered

        // Measurement
        // Use pointer to avoid memory copy
        VectorXd *m;     // Measurement vector
        MatrixXd *V, *H; // Measurement uncertainty matrix V, measurement matrix H

        // State vectors
        VectorXd xp, xf; // xp: predicted, xf: filtered

        // System dynamics
        MatrixXd *F, *Q; // Fi: process matrix, Qi: process uncertainty

        // Uncertainty and residual
        MatrixXd Cp, Cf;         // Variation Matrices. Cp: predicted, Cf: filtered
        VectorXd rp;             // predicted residual
        MatrixXd Rp, Rf, Rp_inv; // residual cov predicted/filtered

        // Chi2
        double chi2_step, chi2_total;
    };
} // namespace kalman

#endif