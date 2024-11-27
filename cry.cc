/*

Copyright (c) 2007-2012, The Regents of the University of California.
Produced at the Lawrence Livermore National Laboratory
UCRL-CODE-227323.
All rights reserved.

For details, see http://nuclear.llnl.gov/simulations
Please also read this http://nuclear.llnl.gov/simulations/additional_bsd.html

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the disclaimer below.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the disclaimer (as noted below) in
the documentation and/or other materials provided with the
distribution.

3. Neither the name of the UC/LLNL nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OF
THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//
//
// write cry particles to text file shower.out
//
#include "CRYGenerator.h"
#include "CRYSetup.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h> // For Ubuntu Linux

// Project includes
#include "libs/cxxopts.hpp"
#include "util.hh"

std::string util::globals::PROJECT_SOURCE_DIR = "";

int main(int argc, const char *argv[])
{
  auto PROJECT_SOURCE_DIR = util::path::getExecutablePath().parent_path().parent_path().string();

  // clang-format off
  cxxopts::Options options("CRY cosmic generator", "CRY cosmic generator with text output");
  options.add_options()
      ("h,help", "Print help")
      ("s,setup", "setup file name", cxxopts::value<std::string>()->default_value("../macros/generators/cry_all.file"))
      ("o,output", "output filename", cxxopts::value<std::string>()->default_value("cry_genparticles.out"))
      ("n,nevents", "Number of events (after energy cut, if energy cut is enabled)", cxxopts::value<int>()->default_value("1000"))
      ("e,ekin", "Kinetic energy cut. Positive value means above, and negative means below. Set to 0 to disable", cxxopts::value<float>()->default_value("0"))
      ("c,contains", "Select only events that contains specific PDG-ids, separated by commas. Eg, --contains=13,-13. All supported: e(+-11), mu(+-13), neutron(2112), proton(2212), gamma(22), pion(+-211,111), kaon(+-321,311,310,130). Set to 0 to disable", cxxopts::value<std::vector<int>>()->default_value("0"));
  auto args = options.parse(argc, argv);
  // clang-format on

  // Show help if the user asks for it
  if (args.count("help"))
  {
    std::cout << options.help() << std::endl;
    return 0; // Exit the program after showing help
  }

  int nEv = args["nevents"].as<int>();
  float ekin_cut = args["ekin"].as<float>();
  int ekin_cut_sign = ekin_cut >= 0 ? 1 : -1;
  auto setup_filename = args["setup"].as<std::string>();
  auto selections = args["contains"].as<std::vector<int>>();
  bool selection_en = selections[0] == 0 ? false : true;
  util::py::print(selections);

  // Apply the setup file
  std::ifstream inputFile;
  inputFile.open(setup_filename, std::ios::in);
  char buffer[1000];
  std::string setupString("");
  while (!inputFile.getline(buffer, 1000).eof())
  {
    setupString.append(buffer);
    setupString.append(" ");
  }
  CRYSetup *setup = new CRYSetup(setupString,
                                 PROJECT_SOURCE_DIR + "/third_party/cry_v1.7/data"); // Use absolute path to CRY data files

  // Make a generator
  CRYGenerator gen(setup);

  // Open output file to write shower particles
  std::ofstream out(args["output"].as<std::string>());
  if (!out || out.bad())
    return 1;
  std::cout << "Output sent to shower.out\n";
  out << "# nEvent nSecondary PDG-id KE[MeV] x[m] y[m] z[m] u v w t[s]\n";

  std::vector<CRYParticle *> *ev = new std::vector<CRYParticle *>;
  for (int i = 0; i < nEv; )
  {
    ev->clear();
      // gen.genEvent(ev);

    if (ekin_cut==0)
      gen.genEvent(ev);
    else
      gen.genEvent(ev, ekin_cut);

    // If cut is enabled, check if the result contains the specified particle
    bool selected = false;
    double tmin = 1e100;
    if (selection_en)
    {
      for (unsigned j = 0; j < ev->size(); j++)
      {
        CRYParticle *p = (*ev)[j];
        // Find the time of the first particle
        if (p->t() < tmin)
          tmin = p->t();
        if ((std::find(selections.begin(), selections.end(), p->PDGid()) != selections.end()) && ((p->ke() * ekin_cut_sign) > ekin_cut))
        {
          selected = true;
          // break;
        }
      }
    }
    else
      selected = true;

    if (!selected)
      continue;
    else
      i++;

    for (unsigned j = 0; j < ev->size(); j++)
    {
      CRYParticle *p = (*ev)[j];
      //      std::cout
      out
          << i                 // event entry number
          << " " << j          // secondary
          << " " << p->PDGid() // particle type
          << " " << p->ke()    // KE
          << " " << p->x()
          << " " << p->y()
          << " " << p->z()
          << " " << p->u()
          << " " << p->v()
          << " " << p->w()
          << " " << p->t() - tmin
          << "\n";
      delete (*ev)[j];
    }
  }

  std::cout << "Run completed.\n";
  std::cout << "Total time simulated: " << gen.timeSimulated() << " seconds\n";
  out << "# Run completed.\n";
  out << "# Total time simulated: " << gen.timeSimulated() << " seconds\n";
  out.close();

  return 0;
}
