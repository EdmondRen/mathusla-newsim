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
/// \file MuSteppingAction.cc
/// \brief Implementation of the MuSteppingAction class

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"


#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"

#include "MuSteppingAction.hh"
#include "MuEventAction.hh"
#include "MuDetectorConstruction.hh"
#include "MuActionInitialization.hh"

// #include "util.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuSteppingAction::MuSteppingAction(
                      const MuDetectorConstruction* detectorConstruction,
                      MuEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuSteppingAction::~MuSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuSteppingAction::UserSteppingAction(const G4Step* step)
{
  // get volume of the current step
  // auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  (void)step;

  const auto step_point = step->GetPreStepPoint();
  // const auto momentum   = G4LorentzVector(step_point->GetTotalEnergy(), step_point->GetMomentum());  
  G4Track* track = (G4Track*)(step->GetTrack());
  G4double Ekin = track->GetKineticEnergy();
  auto _pdg = track->GetParticleDefinition()->GetPDGEncoding();
  auto _trackID = track->GetTrackID();
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();



  // if (MuActionInitialization::DEBUG)
  // {}

  // if ((particleDefinition == G4Electron::Electron() || particleDefinition == G4Positron::Positron()) && Ekin < 10 * MeV) 
  // {
  //   // print("Kill track", _trackID, _pdg, Ekin, step_point->GetTotalEnergy());
  //   track->SetTrackStatus(fStopAndKill);
  // }
  // else if (particleDefinition == G4Gamma::Gamma() && Ekin < 2) 
  // {
  //   // print("Kill track", _trackID, _pdg, Ekin, step_point->GetTotalEnergy());
  //   track->SetTrackStatus(fStopAndKill);
  // }  
  if (particleDefinition == G4Neutron::Neutron() && Ekin < 1 * MeV) 
  {
    // print("Kill track", _trackID, _pdg, Ekin, step_point->GetTotalEnergy());
    track->SetTrackStatus(fStopAndKill);
  }  

  // if (_trackID>10) 
  // {
  //   // print("Kill track", _trackID, _pdg, Ekin, step_point->GetTotalEnergy());
  //   track->SetTrackStatus(fStopAndKill);
  // }  



  // else
    // print(_trackID, _pdg, Ekin, step_point->GetGlobalTime());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
