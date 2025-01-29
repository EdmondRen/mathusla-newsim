// stdlib
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <stdlib.h> // For Ubuntu Linux

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>



// Project includes
#include "libs/cxxopts.hpp"
#include "util.hh"

int main(int argc, const char *argv[])
{
    // Setup argument format
    // clang-format off
    cxxopts::Options options("./tracker", "Reconstruct the track and vertex in the simulation");
    options.add_options()
        ("h,help", "Print help")
        ("filename", "ROOT file to perform reconstruction", cxxopts::value<std::string>())
        ("s,seed", "Seed for random number generator", cxxopts::value<int>()->default_value("-1"))
        ("p,print_progress", "Print progress every `p` events", cxxopts::value<int>()->default_value("1"));
    options.parse_positional({"filename"});
    auto args = options.parse(argc, argv);
    // clang-format on

    // Show help if the user asks for it
    if (args.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0; // Exit the program after showing help
    }

    print("**************************************************************");
    print("   MATHUSLA SIM Tracker, version ", util::VERSION);
    print("      - MATHUSLA Collaboration");
    print("**************************************************************");    


}