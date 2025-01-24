// stdlib
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h> // For Ubuntu Linux


#include "libs/cxxopts.hpp"

// Project includes
#include "util.hh"
#include "generator/parma/parma_util.hh"

// Global variables
std::string util::globals::PROJECT_SOURCE_DIR = "";

int main(int argc, const char *argv[])
{

    // clang-format off
    cxxopts::Options options("CRY cosmic generator", "CRY cosmic generator with text output");
    options.add_options()
        ("h,help", "Print help")
        ("o,output", "output filename", cxxopts::value<std::string>()->default_value("parma_genparticles.out"))
        ("s,setup", "setup file name", cxxopts::value<std::string>()->default_value("../macros/generators/parma_default.conf"))
        ("p,particle", "PDG ID of the particle. Supports n(2112), p(2212), mu+-(-+13), e+-(-+11), gamma(22)", cxxopts::value<int>())
        ("L,ekin_low", "Kinetic energy low end (1 - 2e5)", cxxopts::value<float>())
        ("H,ekin_high", "Kinetic energy high end (1 - 2e5)", cxxopts::value<float>())
        ("n,nevents", "Number of events (after energy cut, if energy cut is enabled)", cxxopts::value<int>()->default_value("1000"));
    auto args = options.parse(argc, argv);
    // clang-format on

    // Show help if the user asks for it
    if (args.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0; // Exit the program after showing help
    }

    int nEv = args["nevents"].as<int>();
    auto setup_filename = args["setup"].as<std::string>();

    // Open output file to write shower particles
    std::ofstream out(args["output"].as<std::string>());
    if (!out || out.bad())
        return 1;
    std::cout << "Output sent to " << args["output"].as<std::string>() << std::endl;
    out << "# nEvent nSecondary PDG-id KE[MeV] x[m] y[m] z[m] u v w t[s]\n";

    // Get the configuration file
    auto parcard = util::io::ParHandler(setup_filename);
    auto config = parcard.GetConfig();



    // Instantiate and config the PARMA generator
    auto generator = PARMA::ParmaGen();

    // Setup the generator with the config map
    generator.configure(config);

    // , or you can manually set all the parameters
    // generator.ip = 1;             // Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
    // generator.iyear = 2019;       // Year
    // generator.imonth = 2;         // Month
    // generator.iday = 1;           // Day
    // generator.glat = 30.5;        // Latitude (deg), -90 =< glat =< 90
    // generator.glong = -76.2;      // Longitude (deg), -180 =< glong =< 180
    // generator.alti = 0.0;         // Altitude (km)
    // generator.g = 0.15;           // Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
    // generator.subboxlength = 100; // Length of the square target area in [cm]
    // generator.emin = 26;        // Minimum energy of particle
    // generator.emax = 1.0e5;       // Maximum energy of particle
    // generator.UpdateParameters();

    // Override some parameter though command line arguments
    if (args.count("ekin_low"))
        generator.emin =  args["ekin_low"].as<float>();        // Minimum energy of particle
    if (args.count("ekin_high"))
        generator.emax =  args["ekin_high"].as<float>();        // Minimum energy of particle  
    if (args.count("particle"))
        generator.ip =  PARMA::pdgid_to_id[args["particle"].as<int>()];        // Minimum energy of particle        
    generator.UpdateParameters();

    for (int i = 0; i < nEv; i++)
    {
        auto p = generator.Generate();
        //      std::cout
        out << i              // event entry number
            << " " << 0       // secondary number, not used here for PARMA
            << " " << p.pdgid // particle type
            << " " << p.ke    // kinetic energy
            << " " << p.x
            << " " << p.y
            << " " << p.z
            << " " << p.u
            << " " << p.v
            << " " << p.w
            << " " << p.t
            << "\n";
    }

    std::cout << "Run completed.\n";
    std::cout << "Flux: " << generator.TotalFlux << " [/cm2/s]\n";
    out << "# Run completed.\n";
    out << "# Flux: " << generator.TotalFlux << " [/cm2/s]\n";
    out.close();

    return 0;
}