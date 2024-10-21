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
/// \file MuDetectorConstruction.cc
/// \brief Implementation of the MuDetectorConstruction class

#include "MuDetectorConstruction.hh"
#include "geometry/_GeoBuilder.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::string MuDetectorConstruction::MessengerDirectory = "/det/";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuDetectorConstruction::MuDetectorConstruction(const std::string &detector_name,
                                               const std::string &export_dir)
    : G4VUserDetectorConstruction(),
      G4UImessenger(MuDetectorConstruction::MessengerDirectory, "Detector constructor"),
      fCheckOverlaps(true)
{
  // Make a map of all available detectors
  _det_map_["uoft1"] = new MuGeoBuilder::Builder("gun", "ParticleGun");
  _det_ = _det_map_["uoft1"];


  // Add messenger commands 
  cmd_select = CreateCommand<G4UIcmdWithAString>("select", "Select Detector.");
  cmd_select->SetParameterName("detector", false);
  cmd_select->SetDefaultValue("uoft1");
  cmd_select->AvailableForStates(G4State_PreInit, G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuDetectorConstruction::~MuDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MuDetectorConstruction::Construct()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuDetectorConstruction::ConstructSDandField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



// -----------------------------------------------------------------------------
// End of G4 requirement
// Start user-define functions
// -----------------------------------------------------------------------------

