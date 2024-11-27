#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <filesystem>

#include "parma_util.hh"

using namespace std;

namespace PARMA
{   
    // Get the path of current folder, and locate PARMA relative to it.
    auto path_current_folder = std::filesystem::path(std::string(__FILE__)).parent_path().string();
    std::string parma_installed_path = path_current_folder + "/../../../external/parma_cpp/";

    const int ParmaGen::IangPart[] = {1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 6};

    // Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
    std::map<int,int> id_to_pdgid = 
       {{0, 2112},
        {1, 2212},
        {29, -13},
        {30, 13},
        {31, 11},
        {32, -11},
        {33, 22}};
    std::map<int,int> pdgid_to_id = 
       {{2112, 0},
        {2212, 1},
        {-13, 29},
        {13, 30},
        {11, 31},
        {-11, 32},
        {22, 33}};

    std::map<int,float> pdgid_to_mass = 
        {{11,0.51099895}, {-11,0.51099895},
         {13,105.6583755},{-13,105.6583755},
         {22,0},
         {2112,940.6},
         {2212,938.27}};



    ParmaGen::ParmaGen()
    {
        UpdateParameters();
        std::cout<< " - PARMA installed at: " << parma_installed_path<<std::endl;
    }

    void ParmaGen::configure(std::map<std::string, double> config)
    {
        std::cout << " ===== PARMA settings ============\n";
        for (auto item:config){
            std::cout<< "  " << item.first<<": " << item.second <<std::endl;
            if (item.first == "particle_pdgid")
                this->ip = pdgid_to_id[static_cast<int>(item.second)];
            else if (item.first == "year")
                this->iyear = static_cast<int>(item.second);
            else if (item.first == "month")
                this->imonth = static_cast<int>(item.second);
            else if (item.first == "day")
                this->iday = static_cast<int>(item.second);
            else if (item.first == "latitude")
                this->glat = item.second;
            else if (item.first == "longitude")
                this->glong = item.second;
            else if (item.first == "altitude")
                this->alti = item.second / 1000; // m to km
            else if (item.first == "water_fraction")
                this->g = item.second;
            else if (item.first == "subboxLength")
                this->subboxlength = item.second * 1000; // m to cm
            else if (item.first == "emin")
                this->emin = item.second;
            else if (item.first == "emax")
                this->emax = item.second;
            else if (item.first == "seed")
                this->seed = item.second;
        }
        std::cout << " ===== END of PARMA settings =====\n";

        UpdateParameters();
    }

    void ParmaGen::UpdateParameters()
    {
        rngdptr = 0;
        setRandomFunction(tmpRandom);

        // calculate parameters
        s = getHPcpp(iyear, imonth, iday); // solar modulation potential
        r = getrcpp(glat, glong);          // Vertical cut-off rigidity (GV)
        d = getdcpp(alti, glat);           // Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

        if (IangPart[ip] == 0)
        {
            cout << "Angular distribution is not available for the particle";
            exit(1);
        }

        if (ip == 0 && emin < 1.0e-8)
        {
            emin = 1.0e-8;
        } // Minimum energy for neutron is 10 meV
        if (ip != 0 && emin < 1.0e-2)
        {
            emin = 1.0e-2;
        } // Minimum energy for other particle is 10 keV

        // Make energy and angle mesh
        double elog = log10(emin);
        double estep = (log10(emax) - log10(emin)) / nebin;
        for (ie = 0; ie <= nebin; ie++)
        {
            ehigh[ie] = pow(10, elog);
            if (ie != 0)
                emid[ie] = sqrt(ehigh[ie] * ehigh[ie - 1]);
            elog = elog + estep;
        }

        // setup for 511 keV annihilation gamma
        if (ip == npart)
        {
            for (ie = 1; ie <= nebin; ie++)
            {
                if (ehigh[ie - 1] < emass && emass <= ehigh[ie])
                {
                    ie511 = ie;
                    flux511bin = getSpecCpp(ip, s, r, d, emid[ie], g) * (ehigh[ie] - ehigh[ie - 1]);
                    annihratio = (flux511bin + get511fluxCpp(s, r, d)) / flux511bin;
                    break;
                }
            }
        }
        else
        {
            ie511 = 0;
        }

        double astep = (amax - amin) / nabin;
        for (ia = 0; ia <= nabin; ia++)
        {
            ahigh[ia] = amin + astep * ia;
            if (ia != 0)
                amid[ia] = (ahigh[ia] + ahigh[ia - 1]) * 0.5;
        }

        // Make probability table (absolute value)
        for (ie = 1; ie <= nebin; ie++)
        {
            for (ia = 1; ia <= nabin; ia++)
            {
                atable[ia][ie] = atable[ia - 1][ie] + getSpecCpp(ip, s, r, d, emid[ie], g) * getSpecAngFinalCpp(IangPart[ip], s, r, d, emid[ie], g, amid[ia]) * (2.0 * acos(-1.0)) * (ahigh[ia] - ahigh[ia - 1]); // angular integrated value
            }
        }
        for (ie = 1; ie <= nebin; ie++)
        {
            if (ip == npart && ie == ie511)
            {
                etable[ie] = etable[ie - 1] + atable[nabin][ie] * (ehigh[ie] - ehigh[ie - 1]) * annihratio; // energy integrated value
            }
            else
            {
                etable[ie] = etable[ie - 1] + atable[nabin][ie] * (ehigh[ie] - ehigh[ie - 1]); // energy integrated value
            }
        }

        TotalFlux = etable[nebin]; // Total Flux (/cm2/s), used for normalization

        // Make probability table (normalized to 1)
        for (ie = 1; ie <= nebin; ie++)
        {
            etable[ie] = etable[ie] / etable[nebin];
            for (ia = 1; ia <= nabin; ia++)
            {
                atable[ia][ie] = atable[ia][ie] / atable[nabin][ie];
            }
        }        
    }

    ParmaParticle ParmaGen::Generate()
    {   
        ParmaParticle particle;

        e = getGenerationCpp(randomFlat(), randomFlat(), &ie, nebin, ehigh, etable); // energy
        if (ip == npart && ie == ie511)
        { // may be annihilation gamma
            if (randomFlat() > 1.0 / annihratio)
            {
                e = emass;
            }
        }
        phi = 2.0 * acos(-1.0) * (randomFlat() - 0.5); // azimuth angle (rad)
        for (ia2 = 0; ia2 <= nabin; ia2++)
        {
            atable2[ia2] = atable[ia2][ie];
        }
        cx = getGenerationCpp(randomFlat(), randomFlat(), &ia, nabin, ahigh, atable2); // z direction, -1.0:upward, 0.0:horizontal, 1.0:downward
        // Sampling position in a circle
        // do
        // {
        //     xd = (randomFlat() - 0.5) * 2.0 * radi;
        //     yd = (randomFlat() - 0.5) * 2.0 * radi;
        // } while (sqrt(xd * xd + yd * yd) > radi);
        // zd = radi;

        // sx = sqrt(1 - cx * cx); // sin(theta)

        // x = xd * cx * cos(phi) - yd * sin(phi) + zd * sx * cos(phi);
        // y = xd * cx * sin(phi) + yd * cos(phi) + zd * sx * sin(phi);
        // z = -xd * sx + zd * cx;
        // u = -sx * cos(phi);
        // v = -sx * sin(phi);
        // w = -cx;

        // sampling position in a square
        x = (randomFlat() - 0.5) * 2.0 * subboxlength;
        y = (randomFlat() - 0.5) * 2.0 * subboxlength;
        z = alti * 1000; // zd in [meter] 
        u = -sx * cos(phi);
        v = -sx * sin(phi);
        w = -cx;

        particle.pdgid = id_to_pdgid[ip];
        particle.ke = e;
        particle.u = u;
        particle.v = v;
        particle.w = w;
        particle.x = x;
        particle.y = y;
        particle.z = z;
        particle.t = 0;

        return particle;
    }

    double ParmaGen::tmpRandom()
    {
        static unsigned long int next = 1;

        next = next * 1103515245 + 123345;
        unsigned temp = (unsigned)(next / 65536) % 32768;

        return (temp + 1.0) / 32769.0;
    }

    // until we get one from the transport coders
    double ParmaGen::randomFlat(double min, double max)
    {
        return min + (max - min) * ((double)rngdptr());
    }


    // ***********************************************************
    double getGenerationCpp(double rand1, double rand2, int *ibin, int nbin, double *high, double *table)
    // ***********************************************************
    {
        double getGeneration = 0.0;

        int i;

        for (i = 1; i <= nbin - 1; i++)
        {
            if (rand1 <= table[i])
            {
                break;
            }
        }
        *ibin = i; // bin ID

        getGeneration = high[*ibin - 1] * rand2 + high[*ibin] * (1.0 - rand2);

        return getGeneration;
    }

} // namespace PARMA
