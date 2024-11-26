#ifndef parma_util_hh
#define parma_util_hh


namespace PARMA
{
    double getHPcpp(int, int, int);
    double getrcpp(double, double);
    double getdcpp(double, double);
    double getSpecCpp(int, double, double, double, double, double);
    double getSpecAngFinalCpp(int, double, double, double, double, double, double);
    double getGenerationCpp(double, double, int *, int, double *, double *);
    double get511fluxCpp(double, double, double);

    // Particle struct
    struct ParmaParticle
    {
        // Necessary parameters
        int pdgid; // PDG identifier
        double px, py, pz;
        double u, v, w;
        double x, y, z, t;

        // Optional parameters
        int index; // A counter, can be initialized to the value you want

        ParmaParticle() = default;

        ParmaParticle(int identifier,
                 double x_momentum,
                 double y_momentum,
                 double z_momentum)
            : pdgid(identifier), px(x_momentum), py(y_momentum), pz(z_momentum), x(0), y(0), z(0), t(0) {}

        ParmaParticle(int identifier,
                 double x_position,
                 double y_position,
                 double z_position,
                 double t_position,
                 double x_momentum,
                 double y_momentum,
                 double z_momentum,
                 int genParticleIndex)
            : pdgid(identifier), px(x_momentum), py(y_momentum), pz(z_momentum),
              x(x_position), y(y_position), z(z_position), t(t_position), index(genParticleIndex) {}
    };    

    class ParmaGen
    {
    public:
        ParmaGen();

        ParmaParticle Generate();

        // Random number generator
        double randomFlat(double min=0., double max=1.);
        static double tmpRandom();
        void setRandomFunction(double (*newFunc)(void)) { rngdptr=newFunc;}



        // Set condition
        int nevent = 1000;    // number of particles to be generated
        int ip = 1;           // Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
        int iyear = 2019;     // Year
        int imonth = 2;       // Month
        int iday = 1;         // Day
        double glat = 30.5;   // Latitude (deg), -90 =< glat =< 90
        double glong = -76.2; // Longitude (deg), -180 =< glong =< 180
        double alti = 0.0;    // Altitude (km)
        double g = 0.15;      // Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
        double radi = 100.0;  // radius of the target area in cm (put your target inside this area)

        // Set energy and angle ranges for generation
        double emin = 1.0e0; // Minimum energy of particle
        double emax = 1.0e5; // Maximum energy of particle
        double amin = -1.0;  // Minimum cosine of particle
        double amax = 1.0;   // Maximum cosine of particle        

        double TotalFlux;

    private:
        double (*rngdptr)(void);        


        const static int npart = 33;
        const static int IangPart[npart + 1];
        const static int nebin = 1000;                          // number of energy mesh (divided by log)
        const static int nabin = 100;                           // number of angle mesh (divided by linear)
        double ehigh[nebin + 1], emid[nebin + 1]; // higher and middle point of energy bin
        double ahigh[nabin + 1], amid[nabin + 1]; // higher and middle point of angular bin
        double etable[nebin + 1];            // probability table (0.0 for 0, 1.0 for nebin)
        double atable[nabin + 1][nebin + 1]; // probability table (0.0 for 0, 1.0 for nabin)
        double atable2[nabin + 1];                       // temporary used dimension for anguluar probability
        double Flux[nabin + 1][nebin + 1];   // Monte Carlo generated flux        

        double emass = 0.51099895e0; // mass of electron in MeV
        double s,r,d;
        int ia, ie, i, ia2, ie511;
        double e, phi, u, v, w, sx, cx, xd, yd, zd, x, y, z, annihratio, flux511bin;

    };

}

#endif