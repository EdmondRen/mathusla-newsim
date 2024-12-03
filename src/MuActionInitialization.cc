//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file MuActionInitialization.cc
/// \brief Implementation of the MuActionInitialization class

#include "MuActionInitialization.hh"
#include "MuPrimaryGeneratorAction.hh"
#include "MuRunAction.hh"
#include "MuEventAction.hh"
#include "MuSteppingAction.hh"
#include "MuDetectorConstruction.hh"

#include "util.hh"

bool MuActionInitialization::ENABLE_DEBUG = false;
bool MuActionInitialization::ENABLE_STEPS = false;


MuActionInitialization::MuActionInitialization(MuDetectorConstruction *detConstruction, cxxopts::ParseResult &uargs)
    : G4VUserActionInitialization(),
      fDetConstruction(detConstruction),
      args(uargs)
{
  // Make new directory for output
  util::io::create_directory(args["output"].as<G4String>());

  // Set static flags
  ENABLE_DEBUG = args["debug"].as<bool>();
  ENABLE_STEPS = args["all_steps"].as<bool>();
}


MuActionInitialization::~MuActionInitialization()
{
}


void MuActionInitialization::BuildForMaster() const
{
  SetUserAction(new MuRunAction(args["output"].as<G4String>(), args["run"].as<G4int>()));
}


void MuActionInitialization::Build() const
{ 
  auto generator_name = args["generator"].as<G4String>();
  auto output_dir = args["output"].as<G4String>();
  auto run_number = args["run"].as<G4int>();
  SetUserAction(new MuPrimaryGeneratorAction(generator_name));
  SetUserAction(new MuRunAction(output_dir,
                                run_number));
  auto eventAction = new MuEventAction;
  SetUserAction(eventAction);
  SetUserAction(new MuSteppingAction(fDetConstruction, eventAction, ENABLE_STEPS));
}

