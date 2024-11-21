/// \file MuPrimaryGeneratorAction.cc
/// \brief Implementation of the MuPrimaryGeneratorAction class


#include "G4RunManager.hh"
#include "Randomize.hh"


// Project include
#include "MuPrimaryGeneratorAction.hh"
#include "util.hh"


MuGenerators::Generator*  _gen_; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuPrimaryGeneratorAction::MuPrimaryGeneratorAction(std::string _gen_default_)
    : G4VUserPrimaryGeneratorAction(), G4UImessenger(MuGenerators::Generator::MessengerDirectory, "Particle Generators.")
{
  // Make a map of all available generators
  _gen_map_["gun"] = new MuGenerators::ParticleGun("gun", "ParticleGun");
  _gen_map_["cry"] = new MuGenerators::MuCRY("cry", "CRY cosmic", util::globals::PROJECT_SOURCE_DIR);
  _gen_ = _gen_map_[_gen_default_];


  // Add messenger commands 
  cmd_select = CreateCommand<G4UIcmdWithAString>("select", "Select Generator.");
  cmd_select->SetParameterName("generator", false);
  cmd_select->SetDefaultValue("gun");
  cmd_select->AvailableForStates(G4State_PreInit, G4State_Idle);
}


MuPrimaryGeneratorAction::~MuPrimaryGeneratorAction()
{
  delete _gen_;
}


void MuPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  // This function is called at the begining of event
  _gen_->GeneratePrimaryVertex(anEvent);
  // std::vector<unsigned long> engienStatus = CLHEP::HepRandom::getTheEngine()->put();
  // util::py::print("RanecuEngine status after gen [address, init_seed, seed[0], seed[1]]", engienStatus);
}


//__Generator Action Messenger Set Value
void MuPrimaryGeneratorAction::SetNewValue(G4UIcommand *command,
                                           G4String value)
{
  if (command == cmd_select){
    _gen_ = _gen_map_[value];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


const MuGenerators::Generator* MuPrimaryGeneratorAction::GetGenerator()
{
  return _gen_;
}

const MuGenerators::ParticleVector MuPrimaryGeneratorAction::GetLastEvent()
{
  return _gen_->genParticles;
}

