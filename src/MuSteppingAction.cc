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
    const MuDetectorConstruction *detectorConstruction,
    MuEventAction *eventAction, bool _ENABLE_STEPS)
    : G4UserSteppingAction(),
      ENABLE_STEPS(_ENABLE_STEPS),fDetConstruction(detectorConstruction), 
      fEventAction(eventAction)
{
  fStepDataStore = new StepDataStore();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuSteppingAction::~MuSteppingAction()
{
  delete fStepDataStore;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuSteppingAction::UserSteppingAction(const G4Step *step)
{
  if (ENABLE_STEPS)
  {
    // fStepDataStore is cleared at the beginning of each event in `MuEventAction`
    const auto step_point = step->GetPreStepPoint();
    const auto track = (G4Track *)(step->GetTrack());
    const auto position4   = G4LorentzVector(step_point->GetGlobalTime(), step_point->GetPosition());
    const auto momentum4   = G4LorentzVector(step_point->GetTotalEnergy(), step_point->GetMomentum());    
    const auto particle = track->GetParticleDefinition();

    fStepDataStore->_step_x.push_back(position4.x());
    fStepDataStore->_step_y.push_back(position4.y());
    fStepDataStore->_step_z.push_back(position4.z());
    fStepDataStore->_step_t.push_back(position4.t());
    fStepDataStore->_step_edep.push_back(step->GetTotalEnergyDeposit());
    fStepDataStore->_step_px.push_back(momentum4.px());
    fStepDataStore->_step_py.push_back(momentum4.py());
    fStepDataStore->_step_pz.push_back(momentum4.pz());
    fStepDataStore->_step_pdg.push_back(particle->GetPDGEncoding());
    fStepDataStore->_step_trackID.push_back(track->GetTrackID());
    fStepDataStore->_step_trackIDparent.push_back(track->GetParentID());
    fStepDataStore->_step_status.push_back(track->GetTrackStatus());
  }

  // get volume of the current step
  // auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
    // Manually kill a track
  // if (particleDefinition == G4Neutron::Neutron() && Ekin < 1 * MeV)
  // {
  //   // print("Kill track", _trackID, _pdg, Ekin, step_point->GetTotalEnergy());
  //   track->SetTrackStatus(fStopAndKill);
  // }

  // if (_trackID>10)
  // {
  //   // print("Kill track", _trackID, _pdg, Ekin, step_point->GetTotalEnergy());
  //   track->SetTrackStatus(fStopAndKill);
  // }
  // print( _trackID, _pdg, Ekin, track->GetGlobalTime(), step_point->GetGlobalTime());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
