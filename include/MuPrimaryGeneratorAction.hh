/// \file MuPrimaryGeneratorAction.hh
/// \brief Definition of the MuPrimaryGeneratorAction class

#ifndef MuPrimaryGeneratorAction_h
#define MuPrimaryGeneratorAction_h 1

// G4 include
#include "G4VUserPrimaryGeneratorAction.hh"
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


// Project include
#include "generator/_Generator.hh"
#include "generator/ParticleGun.hh"
#include "generator/CRY.hh"


class MuPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
public:
  MuPrimaryGeneratorAction();    
  virtual ~MuPrimaryGeneratorAction();

  // Core function to override
  virtual void GeneratePrimaries(G4Event* event);
  
  // set methods
  void SetNewValue(G4UIcommand* command, G4String value);
  static const MuGenerators::Generator* GetGenerator();
  static void SetGenerator(const std::string& generator);  


private:
  // Store all common commands to generator
  G4UIcmdWithAString * cmd_select;

  // The generator in use, _gen_, and a list of generators
  MuGenerators::Generator*  _gen_; 
  std::unordered_map<std::string, MuGenerators::Generator*> _gen_map_;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
