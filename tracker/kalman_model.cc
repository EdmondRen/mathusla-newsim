#include "kalman_model.hh"

namespace Kalman
{

    KalmanTrack4D::KalmanTrack4D(int multiple_scattering,
                                 int ndimMeasure,
                                 int ndimStates) : enMultipleScattering(multiple_scattering),
                                                       Nmeas(ndimMeasure),
                                                       Nstat(ndimStates),
                                                       xf0(Nstat),
                                                       Cf0(Nstat, Nstat),
                                                       Vi(Nmeas-1, Nmeas-1),
                                                       Hi(Nmeas-1, Nstat),
                                                       Fi(Nstat, Nstat), Qi(Nstat, Nstat),
                                                       kf(KF_Forward(Nmeas, Nstat))
    {
        // Measurement matrix is fixed, so we can initialize it here
        Hi(0, 0) = Hi(1, 1) = Hi(2, 2) = 1;
        // Set process matrix to identity
        Fi.setIdentity();
    }

    int KalmanTrack4D::init_state(const Tracker::DigiHit &hit1, const Tracker::DigiHit &hit2)
    {
        Vector3d r1 = hit1.get_vec3();
        Vector3d r2 = hit2.get_vec3();
        Vector3d dr = r2.array() / r1.array();
        float dt = hit2.get_step() - hit1.get_step();
        Vector3d v = dr / dt;

        // Initial State Vector
        this->xf0 << r1(0), r1(1), r1(2), v(0), v(1), v(2);

        // Initial Covariance
        MatrixXd J(6, 8); // Jacobian matrix. Initialized to 0 by default
        MatrixXd err(8, 8);
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
                J(k + 3, j) = -1 / dt;
                J(k + 3, j + 4) = -J(k + 3, j);
                J(k, j + 4) = 1;
                k += 1;
            }
            else
            {
                J(3, j) = v(0) / dt;
                J(4, j) = v(1) / dt;
                J(5, j) = v(2) / dt;
                J(3, j + 4) = -J(3, j);
                J(4, j + 4) = -J(4, j);
                J(5, j + 4) = -J(5, j);
            }
        }
        err.diagonal() << hit1.ex() * hit1.ex(), hit1.ey() * hit1.ey(), hit1.ez() * hit1.ez(), hit1.et() * hit1.et(),
            hit2.ex() * hit2.ex(), hit2.ey() * hit2.ey(), hit2.ez() * hit2.ez(), hit2.et() * hit2.et();
        this->Cf0 = J * err * J.transpose();

        // Record the current step
        this->step_current = hit1.get_step();

        // Initialize the kalman filter instance
        kf.SetInitialState(xf0, Cf0);

        return 0;
    }


    int KalmanTrack4D::new_step(const Tracker::DigiHit &hit)
    {
        this->step_size = hit.get_step() - this->step_current;
        // measurement uncertainty
        this->Vi.diagonal() = hit.get_err3().array().square();
        // dynamics
        this->Fi(3, 0) = this->Fi(4, 1) = this->Fi(5, 2) = this->step_size;
        // multiple scattering
        if (this->enMultipleScattering)
            update_Q(this->step_size);

        // Update kalman filter matrices
        kf.UpdateMatrix(this->Vi, this->Hi, this->Fi, this->Qi);
    }


    int KalmanTrack4D::try_measurement(const Tracker::DigiHit &hit, float sigma_cut, float chi2_cut)
    {
    }

    int KalmanTrack4D::add_measurement(const Tracker::DigiHit &hit)
    {
        kf.Filter(hit.get_vec3());
    }

    Tracker::Track * KalmanTrack4D::run_filter(Tracker::HitList hits)
    {
        init_state(*hits[0], *hits[1]);
        for (int i=1; i<hits.size(); i++)
        {   
            new_step(*hits[i]);
            add_measurement(*hits[i]);
        }

        this->track_recon = new Tracker::Track();

        return this->track_recon;

    }

    int KalmanTrack4D::update_Q(float step, float multiple_scattering_p, float multiple_scattering_length)
    {
    }       

} // namespace Kalman
