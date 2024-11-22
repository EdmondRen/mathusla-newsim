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

// G4 material
#include "G4Material.hh"
#include "G4NistManager.hh"
// G4 geometry
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
// G4 visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
// G4 constants/units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
// G4 manager
#include "G4SDManager.hh"

#include "MuDetectorConstruction.hh"
#include "MuAnalysis.hh"
#include "geometry/_GeoBuilder.hh"
#include "geometry/TestStand_UofT.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::string MuDetectorConstruction::MessengerDirectory = "/det/";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

util::py::Dict *sensitiveDetectorData;
std::string _det_name_;

MuDetectorConstruction::MuDetectorConstruction(const std::string &detector_name,
                                               const std::string &export_dir)
    : G4VUserDetectorConstruction(),
      G4UImessenger(MuDetectorConstruction::MessengerDirectory, "detector construction"),
      fCheckOverlaps(true)
{
  // Make a map of all available detectors
  _det_map_["uoft1"] = new MuGeoBuilder::Uoft1_Builder();
  
  // Make messenger commands
  _ui_select = CreateCommand<G4UIcmdWithAString>("select", "Select Detector to use.");
  _ui_select->SetParameterName("select", false, false);
  _ui_select->SetDefaultValue("uoft1");
  _ui_select->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Set the detector to the default one
  _det_name_ = detector_name;
  _det_ = _det_map_[_det_name_];

  // Private members
  this->_export_dir_ = export_dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuDetectorConstruction::~MuDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MuDetectorConstruction::Construct()
{
  G4cout << "** Constructing detector" << _det_name_ << G4endl;
  _det_->Construct();
  return _det_->worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuDetectorConstruction::ConstructSDandField()
{
  // Assign a sensitive detector to the geometry
  if (_det_name_ == "uoft1")
  {
    auto sensitive_detector = new Analysis::DefaultDetector();
    _det_->ConstructSD(sensitive_detector);
    sensitiveDetectorData = sensitive_detector->GetDataDict();
    G4SDManager::GetSDMpointer()->AddNewDetector(sensitive_detector);
  }
}

util::py::Dict *MuDetectorConstruction::GetSDdata()
{
  return sensitiveDetectorData;
}

const std::string MuDetectorConstruction::GetName()
{
  return _det_name_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuDetectorConstruction::SetNewValue(G4UIcommand *command,
                                         G4String value)
{
  if (command == _ui_select)
  {
    _det_name_ = value;
    _det_ = _det_map_[_det_name_];
    G4cout << "** Selecting detector " << _det_name_ << G4endl;
  }
}
