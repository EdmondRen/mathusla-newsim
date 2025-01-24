
#include "kalman_model.hh"

namespace Kalman
{

    KalmanTrack4D::KalmanTrack4D(int multiple_scattering,
                                 int ndimMeasure,
                                 int ndimStates,
                                 bool debug) : kf(KF_Forward(ndimMeasure - 1, ndimStates)),
                                               DEBUG(debug), enMultipleScattering(multiple_scattering),
                                               Nmeas(ndimMeasure - 1),
                                               Nstat(ndimStates),
                                               xf0(Nstat),
                                               Cf0(Nstat, Nstat),
                                               Vi(Nmeas, Nmeas),
                                               Hi(Nmeas, Nstat),
                                               Fi(Nstat, Nstat), 
                                               Qi(Nstat, Nstat),
                                               track(nullptr)

    {
        Vi.setZero();
        Qi.setZero();
        // Hi.setZero();
        // Measurement matrix is fixed, so we can initialize it here
        Hi.insert(0, 0) = Hi.insert(1, 1) = Hi.insert(2, 2) = 1;
        // Set process matrix to identity
        for (int i = 0; i < Nstat; ++i)
        {
            Fi.insert(i, i) = 1; // Set at (i, i)
        }
    }

    int KalmanTrack4D::init_state(const Tracker::DigiHit &hit1, const Tracker::DigiHit &hit2, bool use_first)
    {
        Vector3d r1 = hit1.get_vec3();
        Vector3d r2 = hit2.get_vec3();
        Vector3d dr = r2.array() - r1.array();
        float dt = hit2.get_step() - hit1.get_step();
        Vector3d v = dr / dt;

        // Initial State Vector
        if (use_first)
            this->xf0 << r1(0), r1(1), r1(2), v(0), v(1), v(2);
        else
            this->xf0 << r2(0), r2(1), r2(2), v(0), v(1), v(2);

        // Initial Covariance
        MatrixXdSp J(6, 8); // Jacobian matrix. Initialized to 0 by default
        MatrixXdSp err(8, 8);
        // Set element of the Jacobian. Example:
        // [ 0       , 0           , 0       , 0       , 1       , 0             , 0     , 0     ],
        // [ 0       , 0           , 0       , 0       , 0       , 0             , 1     , 0     ],
        // [ 0       , 0           , 0       , 0       , 0       , 0             , 0     , 1     ],
        // [- 1 / dy, dx / (dy*dy) , 0       , 0       , 1 / dy  , - dx / (dy*dy), 0     , 0     ],
        // [0       , dz / (dy*dy) , - 1 / dy, 0       , 0       , - dz / (dy*dy), 1 / dy, 0     ],
        // [0       , dt / (dy*dy) , 0       , - 1 / dy, 0       , - dt / (dy*dy), 0     , 1 / dy]]
        for (int j = 0, k = 0; j < 4; ++j)
        {
            if (j != hit1.param_ind)
            {
                J.insert(k + 3, j) = -1 / dt;
                J.insert(k + 3, j + 4) = -J.coeff(k + 3, j);
                if (use_first)
                    J.insert(k, j) = 1;
                else
                    J.insert(k, j + 4) = 1;
                k += 1;
            }
            else
            {
                J.insert(3, j) = v(0) / dt;
                J.insert(4, j) = v(1) / dt;
                J.insert(5, j) = v(2) / dt;
                J.insert(3, j + 4) = -J.coeff(3, j);
                J.insert(4, j + 4) = -J.coeff(4, j);
                J.insert(5, j + 4) = -J.coeff(5, j);
            }
        }
        err.diagonal() << hit1.ex() * hit1.ex(), hit1.ey() * hit1.ey(), hit1.ez() * hit1.ez(), hit1.et() * hit1.et(),
            hit2.ex() * hit2.ex(), hit2.ey() * hit2.ey(), hit2.ez() * hit2.ez(), hit2.et() * hit2.et();
        this->Cf0 = J * err * J.transpose();

        // Record the position of current step
        this->step_current = use_first ? hit1.get_step() : hit2.get_step();

        // Initialize the kalman filter instance
        kf.SetInitialState(xf0, Cf0);

        if (DEBUG)
            std::cout << "Tracker: ->Kalman filter initialized with\n    ->state vector:\n"
                      << xf0.transpose() << "\n    ->covariance:\n"
                      << Cf0 << std::endl;

        return 0;
    }

    int KalmanTrack4D::new_step(const Tracker::DigiHit &hit)
    {
        this->step_size = hit.get_step() - this->step_current;
        this->step_current = hit.get_step();
        // measurement uncertainty
        this->Vi.diagonal() = hit.get_err3().array().square();
        // dynamics
        this->Fi.coeffRef(0, 3) = this->Fi.coeffRef(1, 4) = this->Fi.coeffRef(2, 5) = this->step_size;

        // multiple scattering
        if (this->enMultipleScattering)
            update_Q(this->step_size);

        // Update kalman filter matrices
        kf.UpdateMatrix(this->Vi, this->Hi, this->Fi, this->Qi);

        if (DEBUG)
        {
            std::cout << "Tracker: -> New step with step size " << step_size << std::endl;
            std::cout << "V (measurement covariance):\n"
                      << Vi << std::endl;
            std::cout << "H (measurement matrix):\n"
                      << Hi.toDense() << std::endl;
            std::cout << "F (state propogation matrix):\n"
                      << Fi.toDense() << std::endl;
        }

        return 0;
    }

    int KalmanTrack4D::try_measurement(const Tracker::DigiHit &hit, float sigma_cut, float chi2_cut)
    {
        (void)hit;
        (void)sigma_cut;
        (void)chi2_cut;
        return 0;
    }

    int KalmanTrack4D::add_measurement(const Tracker::DigiHit &hit)
    {

        kf.Filter(hit.get_vec3());
        if (DEBUG)
        {
            std::cout << "Tracker: -> Add measurement: \n " << hit.vec4.transpose() << std::endl;
            std::cout << "Tracker: -> New filtered state: \n " << kf.GetState().transpose() << std::endl;
            std::cout << "         -> chi2 contribution: " << kf.GetChi2Step() << std::endl;
        }
        return 0;
    }

    std::unique_ptr<Tracker::Track> KalmanTrack4D::run_filter(const std::vector<std::unique_ptr<Tracker::DigiHit>> &hits)
    {
        bool use_first_hit = false;
        init_state(*hits[0], *hits[1], use_first_hit);
        for (size_t i = 2; i < hits.size(); i++)
        {
            new_step(*hits[i]);
            add_measurement(*hits[i]);
        }

        // Insert the independent variable back to the parameter list and covariance matrix
        // After the insertion, there will always be 8 parameters in sequence: {x0, y0, z0, t0, Ax, Ay, Az, At}
        auto indep_param_idx = hits.back()->param_ind;
        auto params = Tracker::Helper::insertVector(kf.GetState(), indep_param_idx, hits.back()->get_step());
        auto cov = Tracker::Helper::insertRowAndColumnOfZeros(kf.GetCov(), indep_param_idx);
        cov = Tracker::Helper::insertRowAndColumnOfZeros(cov, indep_param_idx + 4);
        cov(indep_param_idx, indep_param_idx) = std::pow(hits.back()->get_steperr(), 2);

        // Make the track
        track = std::make_unique<Tracker::Track>();
        track->params = params;
        track->param_ind = indep_param_idx;
        track->cov = cov;
        track->chi2 = kf.GetChi2();

        if (DEBUG)
        {
            std::cout << "Tracker: -> track finished, independent variable added to the parameter list." << std::endl;
            std::cout << "Final filtered state vec:\n"
                      << track->params << std::endl;
            std::cout << "Final filtered state cov:\n"
                      << track->cov << std::endl;
            std::cout << "Final filtered chi2:\n"
                      << track->chi2 << std::endl;
        }

        return std::move(track);
    }

    int KalmanTrack4D::update_Q(float step, float multiple_scattering_p, float multiple_scattering_length)
    {
        (void)step;
        (void)multiple_scattering_p;
        (void)multiple_scattering_length;
        return 0;
    }

} // namespace Kalman
