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
/// \file MuRunAction.cc
/// \brief Implementation of the MuRunAction class
#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <thread>
#include <fstream>
#include <ostream>

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// #include "TFile.h"
// #include "TNamed.h"
// #include "TTree.h"
// #include "TChain.h"

#include "MuRunAction.hh"
#include "MuAnalysis.hh"
#include "MuDetectorConstruction.hh"
#include "MuPrimaryGeneratorAction.hh"

//----------------------------------------------------------------------------------------------
// Helper functions

// template <class... Args>
// void _write_setting(const std::string &name, Args &&...args)
// {
//   auto content = util::py::fprint(args...);
//   TNamed setting(name.c_str(), content.c_str());
//   setting.Write();
// }

std::string _get_datetime()
{
  // Get the current time as a time_point
  auto now = std::chrono::system_clock::now();

  // Convert to a time_t for formatting
  std::time_t now_c = std::chrono::system_clock::to_time_t(now);

  // Format the time as a string
  std::ostringstream oss;
  oss << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
  return oss.str();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuRunAction::MuRunAction(std::string OutDir, int RunNumber)
    : G4UserRunAction(), output_dir(OutDir), run_number(RunNumber)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in MuAnalysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuRunAction::~MuRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuRunAction::BeginOfRunAction(const G4Run * /*run*/)
{
  // inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  const auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  // Analysis::Open(fileName);
  this->fileName_output = output_dir + "/run_" + std::to_string(run_number) + ".root";
  analysisManager->OpenFile(this->fileName_output);

  // Creat tuple based on the sensitive detector of a geometry
  // auto detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  // auto thisDetectorConstruction = dynamic_cast<const MuDetectorConstruction*>(detectorConstruction);
  auto sensitiveDetectorData = MuDetectorConstruction::GetSDdata();
  if (sensitiveDetectorData)
    Analysis::MuCreateNTuple(*sensitiveDetectorData, "raw"); // Tree name is "raw"


  // Create another ntuple for metadata
  this->tupleID_metadata = analysisManager->CreateNtuple("metadata", "Simulation metadata");
  analysisManager->CreateNtupleSColumn(tupleID_metadata, "SimulationName"); // 0
  analysisManager->CreateNtupleSColumn(tupleID_metadata, "Geometry"); // 1
  analysisManager->CreateNtupleSColumn(tupleID_metadata, "Generator"); // 2
  analysisManager->CreateNtupleSColumn(tupleID_metadata, "Time"); // 3
  analysisManager->FinishNtuple(tupleID_metadata);

  // Write metadata
  // analysisManager->FillNtupleSColumn(tupleID_metadata, 0, util::py::fprint("Mathusla simulation"));
  // analysisManager->FillNtupleSColumn(tupleID_metadata, 1, util::py::fprint(MuDetectorConstruction::GetName()));
  // analysisManager->FillNtupleSColumn(tupleID_metadata, 2, util::py::fprint(MuPrimaryGeneratorAction::GetName()));
  // analysisManager->FillNtupleSColumn(tupleID_metadata, 3, util::py::fprint(_get_datetime()));
  
  // analysisManager->AddNtupleRow(tupleID_metadata);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MuRunAction::EndOfRunAction(const G4Run * /*run*/)
{

  auto analysisManager = G4AnalysisManager::Instance();

  // Write metadata
  analysisManager->FillNtupleSColumn(tupleID_metadata, 0, util::py::fprint("Mathusla simulation"));
  analysisManager->FillNtupleSColumn(tupleID_metadata, 1, util::py::fprint(MuDetectorConstruction::GetName()));
  analysisManager->FillNtupleSColumn(tupleID_metadata, 2, util::py::fprint(MuPrimaryGeneratorAction::GetName()));
  analysisManager->FillNtupleSColumn(tupleID_metadata, 3, util::py::fprint(_get_datetime()));
  
  analysisManager->AddNtupleRow(tupleID_metadata);  

  // Save the ROOT file
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
