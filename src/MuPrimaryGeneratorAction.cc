/// \file MuPrimaryGeneratorAction.cc
/// \brief Implementation of the MuPrimaryGeneratorAction class


#include "G4RunManager.hh"
#include "Randomize.hh"


// Project include
#include "MuPrimaryGeneratorAction.hh"
#include "util.hh"
// Global variables
#include "util_globals.hh"

#include "generator/ParticleGun.hh"
#include "generator/CRY.hh"
#include "generator/PARMA.hh"
#include "generator/Recreate.hh"


MuGenerators::Generator*  _gen_; 
std::string _gen_name_;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MuPrimaryGeneratorAction::MuPrimaryGeneratorAction(std::string _gen_default_)
    : G4VUserPrimaryGeneratorAction(), G4UImessenger(MuGenerators::Generator::MessengerDirectory, "Particle Generators.")
{
  // Make a map of all available generators
  _gen_map_["gun"] = new MuGenerators::ParticleGun("gun", "ParticleGun");
  _gen_map_["cry"] = new MuGenerators::MuCRY("cry", "CRY cosmic generator", util::globals::PROJECT_SOURCE_DIR);
  _gen_map_["parma"] = new MuGenerators::MuPARMA("parma", "PARMA cosmic generator", util::globals::PROJECT_SOURCE_DIR);
  _gen_map_["recreate"] = new MuGenerators::MuRecreate("recreate", "Recreate generator, rerun previous events", util::globals::PROJECT_SOURCE_DIR);
  _gen_name_ = _gen_default_;
  _gen_ = _gen_map_[_gen_name_];


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
    _gen_name_ = value;
    _gen_ = _gen_map_[value];
    // print("_gen_name_ is: ", _gen_name_);
    // print("MuPrimaryGeneratorAction::GetName() is: ", MuPrimaryGeneratorAction::GetName());

    // if (_gen_name_=="recreate")
    // {
    // G4UImanager* UImanager = G4UImanager::GetUIpointer();
    // if(UImanager->GetUIpointer()) {
    //     UImanager->ApplyCommand("/run/beamOn 1000");
    // }
    // }

  }
}

G4String MuPrimaryGeneratorAction::GetCurrentValue(G4UIcommand* command)
{
  if (command == cmd_select){
    return _gen_name_;
  }
  
  return "";
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


const MuGenerators::Generator* MuPrimaryGeneratorAction::GetGenerator()
{
  return _gen_;
}

const std::string MuPrimaryGeneratorAction::GetName()
{
  return _gen_name_;
}

const MuGenerators::ParticleVector MuPrimaryGeneratorAction::GetLastEvent()
{
  return _gen_->genParticles;
}

