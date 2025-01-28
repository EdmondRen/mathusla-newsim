
#include "kalman_model.hh"
#include "util.hh"

namespace Kalman
{
    // ------------------------------------------------------------------------------------------------------
    // KalmanTrack4D

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
            if (j != hit1.iv_index)
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

    float KalmanTrack4D::try_measurement(const Tracker::DigiHit &hit, float sigma_cut, float chi2_cut)
    {
        // First, check if the hit is within the predicted range
        Vector3d mi = hit.get_vec3();
        // Vector3d mp = kf.GetPredict().array();
        // Vector3d mperr = kf.GetPredictErr().array();
        // if (std::abs(mi.x()-mp(0)) < sigma_cut*mperr(0))

        // Then, check if it have small enough chi2
        auto chi2_temp = kf.TryMeasurement(mi, sigma_cut);
        if (chi2_temp < chi2_cut)
            return chi2_temp;
        else
            return -1;
    }

    int KalmanTrack4D::add_measurement(const Tracker::DigiHit &hit)
    {

        this->step_current = hit.get_step();
        kf.Filter(hit.get_vec3());
        if (DEBUG)
        {
            std::cout << "Tracker: -> Add measurement: \n " << hit.vec4.transpose() << std::endl;
            std::cout << "Tracker: -> New filtered state: \n " << kf.GetState().transpose() << std::endl;
            std::cout << "         -> chi2 contribution: " << kf.GetChi2Step() << std::endl;
        }
        return 0;
    }

    std::unique_ptr<Tracker::Track> KalmanTrack4D::run_filter(const std::vector<Tracker::DigiHit *> &hits)
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
        auto indep_variable_idx = hits.back()->iv_index;
        auto params_full = Tracker::Helper::insertVector(kf.GetState(), indep_variable_idx, hits.back()->get_step());
        params_full = Tracker::Helper::insertVector(params_full, indep_variable_idx + 4, 1);
        auto cov_full = Tracker::Helper::insertRowAndColumnOfZeros(kf.GetCov(), indep_variable_idx);
        cov_full = Tracker::Helper::insertRowAndColumnOfZeros(cov_full, indep_variable_idx + 4);
        cov_full(indep_variable_idx, indep_variable_idx) = std::pow(hits.back()->get_steperr(), 2);
        cov_full(indep_variable_idx + 4, indep_variable_idx + 4) = 0.0001;

        // Make the track
        track = std::make_unique<Tracker::Track>();
        track->params = kf.GetState();
        track->cov = kf.GetCov();
        track->chi2 = kf.GetChi2();
        track->iv_index = indep_variable_idx;
        track->iv_value = hits.back()->get_step();
        track->iv_error = hits.back()->get_steperr();
        track->params_full = params_full;
        track->cov_full = cov_full;

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

    // ------------------------------------------------------------------------------------------------------
    // LSVertex4DFitter

    std::vector<Tracker::Track *> LSVertex4DFitter::tracks = {};

    void LSVertex4DFitter::cost(int &npar, double *gin, double &f, double *par, int iflag)
    {
        (void)npar;
        (void)gin;
        (void)iflag;

        double chi2_total = 0;
        Vector4d vertex(par[0], par[1], par[2], par[3]); // Fetch the parameters {x,y,z,t}

        for (auto track : LSVertex4DFitter::tracks)
        {
            auto dist_chi2 = track->get_same_time_dist_chi2(vertex);
            chi2_total += dist_chi2.second;
        }

        f = chi2_total; // There is no return value. Return is passed by argument "f"
    }

    bool LSVertex4DFitter::fit(std::vector<Tracker::Track *> _tracks, std::vector<double> arg_guess, double tolerance, double maxcalls)
    {
        LSVertex4DFitter::tracks = _tracks;

        std::vector<double> guess(4);
        if (arg_guess.size() != 0)
            guess = arg_guess;
        else
        {
            auto pair = Tracker::Track::get_closest_midpoint(*tracks[0], *tracks[1]);
            std::copy(pair.second.data(), pair.second.data() + pair.second.size(), guess.begin());
        }

        // Configure the parameters {index, name, initial_value, initial_step, bound_low, bound_high}
        int ierflg = 0; // output error flag
        double first_step_size_x = 1;
        double first_step_size_t = 0.01;
        minimizer.mnparm(0, "x", guess[0], first_step_size_x, 0, 0, ierflg);
        minimizer.mnparm(1, "y", guess[1], first_step_size_x, 0, 0, ierflg);
        minimizer.mnparm(2, "z", guess[2], first_step_size_x, 0, 0, ierflg);
        minimizer.mnparm(3, "t", guess[3], first_step_size_t, 0, 0, ierflg);

        // Run the minimizer
        double arglist[2] = {maxcalls, tolerance};
        minimizer.mnexcm("MIGRAD", arglist, 2, ierflg);

        // Get fit results
        for (int ii = 0; ii < npar; ii++)
        {
            minimizer.GetParameter(ii, parameters[ii], parameter_errors[ii]);
            params(ii) = parameters[ii];
        }
        // get covariance matrix
        minimizer.mnemat(&cov_matrix[0][0], npar);
        // Copy the data from the C-style array to Eigen matrix
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
                this->cov(i, j) = cov_matrix[i][j];
        }

        // Get minimizer status
        double fedm = 0.0;
        double errdef = 0.0;
        int npari = 0;
        int nparx = 0;
        // status:
        // 0 - no covariance
        // 1 - not accurate, diagonal approximation
        // 2 - forced positive definite
        // 3 - full accurate matrix--succesful convergence
        minimizer.mnstat(this->chi2, fedm, errdef, npari, nparx, this->status);

        return (this->status >= 2) ? true : false;
    }

    // ------------------------------------------------------------------------------------------------------
    // KalmanVertex4D

    KalmanVertex4D::KalmanVertex4D(bool debug) : kf(KF_Forward(4, 4)),
                                                 DEBUG(debug), enMultipleScattering(true),
                                                 Nmeas(4),
                                                 Nstat(4),
                                                 xf0(Nstat),
                                                 Cf0(Nstat, Nstat),
                                                 Vi(Nmeas, Nmeas),
                                                 Hi(Nmeas, Nstat),
                                                 Fi(Nstat, Nstat),
                                                 Qi(Nstat, Nstat),
                                                 vertex(nullptr)

    {
        Vi.setZero();
        Qi.setZero();
        // Hi.setZero();
        // Measurement matrix is just identity matrix, so we can initialize it here
        Hi.setIdentity();
        // Process matrix to identity
        Fi.setIdentity();
    }

    int KalmanVertex4D::init_state(VertexSeed &seed)
    {
        // Initial State Vector
        this->xf0 = seed.vertex_fit;

        // Initial Covariance
        this->Cf0 = seed.cov;

        // Initialize the kalman filter instance
        kf.SetInitialState(xf0, Cf0);

        if (DEBUG)
            std::cout << "Tracker: ->Kalman filter initialized with\n    ->state vector:\n"
                      << xf0.transpose() << "\n    ->covariance:\n"
                      << Cf0 << std::endl;

        return 0;
    }

    float KalmanVertex4D::try_measurement(const Tracker::Track &track, float distance_cut, float chi2_cut)
    {
        // measurement uncertainty
        auto meas_quality = track.get_closest_point_and_cov(this->kf.GetState());
        this->Vi = meas_quality.second;
        this->mi = meas_quality.first;
        auto distance_to_vertex = (this->mi-this->kf.GetState()).segment(0,3).norm();

        // dynamics: Do nothing, Fi matrix is always identity.

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

        // Then, check if it have small enough chi2
        auto chi2_temp = kf.TryMeasurement( this->mi, -1);
        if (chi2_temp < chi2_cut && distance_to_vertex < distance_cut)
            return chi2_temp;
        else
            return -1;
    }

    int KalmanVertex4D::add_measurement()
    {

        kf.Filter(this->mi);
        if (DEBUG)
        {
            std::cout << "Tracker: -> Add measurement: \n " << this->mi.transpose() << std::endl;
            std::cout << "Tracker: -> New filtered state: \n " << kf.GetState().transpose() << std::endl;
            std::cout << "         -> chi2 contribution: " << kf.GetChi2Step() << std::endl;
        }
        return 0;
    }

    std::unique_ptr<Tracker::Vertex> KalmanVertex4D::run_filter(const std::vector<Tracker::Track *> &tracks, VertexSeed *seed)
    {   
        std::unique_ptr<VertexSeed> seed_temp;
        if (seed==nullptr)
        {
            seed_temp = std::make_unique<VertexSeed>(tracks[0], tracks[1]);
            seed_temp->GetScore();
            seed = seed_temp.get();
        }

        init_state(*seed);
        for (size_t i = 2; i < tracks.size(); i++)
        {
            try_measurement(*tracks[i], 4000, 10);
            add_measurement();
        }

        // Make the track
        vertex = std::make_unique<Tracker::Vertex>();
        vertex->params = kf.GetState();
        vertex->cov = kf.GetCov();
        vertex->chi2 = kf.GetChi2();

        if (DEBUG)
        {
            std::cout << "Tracker: -> Vertex finished, independent variable added to the parameter list." << std::endl;
            std::cout << "Final filtered state vec:\n"
                      << vertex->params << std::endl;
            std::cout << "Final filtered state cov:\n"
                      << vertex->cov << std::endl;
            std::cout << "Final filtered chi2:\n"
                      << vertex->chi2 << std::endl;
        }

        return std::move(vertex);
    }

    // int KalmanVertex4D::update_Q(float step, float multiple_scattering_p, float multiple_scattering_length)
    // {
    //     (void)step;
    //     (void)multiple_scattering_p;
    //     (void)multiple_scattering_length;
    //     return 0;
    // }

} // namespace Kalman
