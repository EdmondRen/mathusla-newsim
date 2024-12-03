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
/// \file MuEventAction.cc
/// \brief Implementation of the MuEventAction class

#include "MuEventAction.hh"
#include "MuRunAction.hh"
#include "MuAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

// # Project include 
#include "MuSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Event information class
// Constructor
MyEventInformation::MyEventInformation(unsigned long seed,
                                       unsigned long seed0,
                                       unsigned long seed1) : fSeed_init(seed),
                                                              fSeed_0(seed0),
                                                              fSeed_1(seed1) {}

// Getter for the seed
std::vector<unsigned long> MyEventInformation::GetInfo() const
{
  std::vector<unsigned long> v;
  v.push_back(fSeed_init);
  v.push_back(fSeed_0);
  v.push_back(fSeed_1);
  return v;
}

// Override Print method
void MyEventInformation::Print() const
{
  std::cout << "Seed_init: " << fSeed_init << std::endl;
  std::cout << "Seed 0: " << fSeed_0 << std::endl;
  std::cout << "Seed 1: " << fSeed_1 << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// EventAction class

MuEventAction::MuEventAction()
    : G4UserEventAction()
{
}

MuEventAction::~MuEventAction()
{
}

void MuEventAction::BeginOfEventAction(const G4Event *event)
{
  // initialisation per event

  // Get the data related to random engine status
  // G4Random::showEngineStatus(); // Directly show it at stdout

  // 1. This only works for RanecuEngine engine!
  // 2. RanecuEngine engine status is described by three unsigned ints
  //    {init_seed, seed[0], seed[1]}
  //    Any other engines are more complicated than this.
  // 3. The recorded status is AFTER the generator finished!
  //    So that it can be used to restore the run with recorded generated particles.
  std::vector<unsigned long> engienStatus = CLHEP::HepRandom::getTheEngine()->put();
  // util::py::print(" Begin of event seed RanecuEngine status [address, init_seed, seed[0], seed[1]]", engienStatus);

  auto *eventInfo = new MyEventInformation(engienStatus[1], engienStatus[2], engienStatus[3]);
  const_cast<G4Event *>(event)->SetUserInformation(eventInfo);

  // Clear the step data store
  auto userStepAction = dynamic_cast<const MuSteppingAction *>(G4RunManager::GetRunManager()->GetUserSteppingAction());
  userStepAction->fStepDataStore->clear();

}

void MuEventAction::EndOfEventAction(const G4Event *event)
{
  // get analysis manager
  // auto analysisManager = G4AnalysisManager::Instance();

  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ((printModulo > 0) && (eventID % printModulo == 0))
  {
    G4cout << "---> End of event: " << eventID << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
