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
/// \file MuEventAction.hh
/// \brief Definition of the MuEventAction class

#ifndef MuEventAction_h
#define MuEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4VUserEventInformation.hh"

#include "globals.hh"


// Custom event information class
// Use this to save per-event information
// * Random number seed
class MyEventInformation : public G4VUserEventInformation {
public:
    // Constructor
    explicit MyEventInformation(long seed);

    // Getter for the seed
    long GetSeed() const;

    // Overriding the pure virtual function from G4VUserEventInformation
    virtual void Print() const override;

    // Virtual destructor
    virtual ~MyEventInformation() override = default;

private:
    long fSeed; // Store the random seed
};



/// Event action class


class MuEventAction : public G4UserEventAction
{
  public:
    MuEventAction();
    virtual ~MuEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddAbs(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);
    
  private:
    int counter;
};

// inline functions

                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
