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
/// \file MuDetectorConstruction.hh
/// \brief Definition of the MuDetectorConstruction class

#ifndef MuDetectorConstruction_h
#define MuDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <G4UIcommand.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UImessenger.hh>
#include <G4UImanager.hh>
#include "G4GenericMessenger.hh"

#include "geometry/_GeoBuilder.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class MuDetectorConstruction : public G4VUserDetectorConstruction, public G4UImessenger
{
  public:
    MuDetectorConstruction(const std::string& detector_name,
          const std::string& export_dir);
    virtual ~MuDetectorConstruction();

  public:
    // Core function to override
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    
    // Handel messenger commands

    // get methods
     
    static const std::string MessengerDirectory;

  private:
    // methods
  
    // data members
    // The detector in use, _det_, with name _det_name_
    // and a list of detectors
    MuGeoBuilder::Builder*  _det_; 
    std::string _det_name_;
    std::unordered_map<std::string, MuGeoBuilder::Builder*> _det_map_;    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    // Store messenger related variables
    void SetNewValue(G4UIcommand *command,
                                G4String value) override;    
    // Messenger commands
    G4UIcmdWithAString *_ui_select;
};

// inline functions


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif