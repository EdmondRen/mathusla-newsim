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





#include "MuDetectorConstruction.hh"

#include "geometry/_GeoBuilder.hh"
#include "geometry/TestStand_UofT.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::string MuDetectorConstruction::MessengerDirectory = "/det/";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuDetectorConstruction::MuDetectorConstruction(const std::string &detector_name,
                                               const std::string &export_dir)
    : G4VUserDetectorConstruction(),
      fCheckOverlaps(true)
{
  // Make a map of all available detectors
  _det_map_["uoft1"] = new MuGeoBuilder::Uoft1_Builder();

  // Add messenger commands 
  fMessenger  = new G4GenericMessenger(this, MuDetectorConstruction::MessengerDirectory, "Detector constructor"),
  fMessenger->DeclareProperty("select", mes_det_selected, "Select Detector to use.");
  mes_det_selected = detector_name; // Default value

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuDetectorConstruction::~MuDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MuDetectorConstruction::Construct()
{
  G4cout<<"Construct start -------"<<mes_det_selected<< G4endl;
  _det_map_[mes_det_selected]->Construct();
  return _det_map_[mes_det_selected]->worldPV;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuDetectorConstruction::ConstructSDandField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

