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
/// \file MuSteppingAction.hh
/// \brief Definition of the MuSteppingAction class

#ifndef MuSteppingAction_h
#define MuSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class MuDetectorConstruction;
class MuEventAction;

class StepDataStore
{
public:
  std::vector<float> _step_x;
  std::vector<float> _step_y;
  std::vector<float> _step_z;
  std::vector<float> _step_t;
  std::vector<float> _step_edep;
  std::vector<float> _step_px;
  std::vector<float> _step_py;
  std::vector<float> _step_pz;
  std::vector<uint32_t> _step_trackID;
  std::vector<uint32_t> _step_trackIDparent;
  std::vector<uint32_t> _step_pdg;
  std::vector<uint32_t> _step_status;

  void clear()
  {
    _step_x.clear();
    _step_y.clear();
    _step_z.clear();
    _step_t.clear();
    _step_edep.clear();
    _step_px.clear();
    _step_py.clear();
    _step_pz.clear();
    _step_pdg.clear();
    _step_trackID.clear();
    _step_trackIDparent.clear();
    _step_status.clear();
  }
};

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track
/// lengths of charged particles in Absober and Gap layers and
/// updated in MuEventAction.

class MuSteppingAction : public G4UserSteppingAction
{
public:
  MuSteppingAction(const MuDetectorConstruction *detectorConstruction,
                   MuEventAction *eventAction, bool ENABLE_STEPS);
  virtual ~MuSteppingAction();

  virtual void UserSteppingAction(const G4Step *step);

  bool ENABLE_STEPS;
  StepDataStore *fStepDataStore;

private:
  const MuDetectorConstruction *fDetConstruction;
  MuEventAction *fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
