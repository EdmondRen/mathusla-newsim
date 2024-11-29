/// \file simulation.cc
/// \brief Main program of the mathusla simulation

// Geant4 include
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4RunManagerFactory.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4StepLimiterPhysics.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh"

// Project include
#include "MuDetectorConstruction.hh"
#include "MuActionInitialization.hh"
#include "MuPrimaryGeneratorAction.hh"
#include "util.hh"
#include "libs/cxxopts.hpp"

#include <sys/stat.h>

// Global variables
std::string util::globals::PROJECT_SOURCE_DIR = "";

int main(int argc, char **argv)
{

  // Get the source file directory
  auto PROJECT_SOURCE_DIR = util::path::getExecutablePath().parent_path().parent_path().string();
  G4cout << PROJECT_SOURCE_DIR << G4endl;
  util::globals::PROJECT_SOURCE_DIR = PROJECT_SOURCE_DIR;

  // ------------------------------
  // Evaluate arguments
  // clang-format off
  cxxopts::Options options("Geant4 simulation", "Geant4 simulation for MATHUSLA");
  options.add_options()
      ("h,help", "Print help")
      ("D,debug", "Enable debugging") // a bool parameter
      ("m,macro", "Macro name, followed by possible parameters for macro seperated by commas. For example, -m=run1.mac,Ek,10,theta,20", cxxopts::value<std::vector<std::string>>())
      ("g,generator", "Generator, one of <gun/range/parma/cry/filereader>", cxxopts::value<G4String>()->default_value("gun"))
      ("d,detector", "Detector, one of <math40/uoft1>", cxxopts::value<G4String>()->default_value("uoft1"))
      ("e,export", "Export directory for geometry and other information. Default is ./export/", cxxopts::value<G4String>()->default_value("export"))
      ("o,output", "Output directory. Default is ./data/", cxxopts::value<G4String>()->default_value("data"))
      ("r,run", "Run number", cxxopts::value<G4int>()->default_value("0"))
      ("s,seed", "Seed of random number generator, a positive integer. Default to -1 for random seed. Events in the same run share the same seed.", cxxopts::value<G4int>()->default_value("-1"))
      ("S,session", "Session name", cxxopts::value<G4String>()->default_value("MathuslaSim"))
      ("t,threads", "[NOT IMPLEMENTED] Number of threads. Only works with single (1) thread at this moment.", cxxopts::value<G4int>()->default_value("1"));
  auto args = options.parse(argc, argv);
  // Show help if the user asks for it
  if (args.count("help"))
  {
    std::cout << options.help() << std::endl;
    return 0; // Exit the program after showing help
  }
  // clang-format on
  // ------------------------------

  // Read and process arguments
  // * macro: need to separate the macro_path and forwarding arguments
  //   The first element of macro_commands is the macro path, the rest are key,value pairs.
  //   For example, a macro_commands could be {run1.mac, Ek, 10, theta, 20}
  std::vector<std::string> macro_commands;
  std::string macro_path;
  macro_commands = args["macro"].count() ? args["macro"].as<std::vector<std::string>>() : macro_commands;
  if (macro_commands.size())
    macro_path = macro_commands[0];

  // Setup random number generator
  G4int run_number = args["run"].as<G4int>();
  (void)run_number;
  G4int run_seed = args["seed"].as<G4int>();
  // Optionally: choose a different Random engine (default is MixMaxRng)
  // Tom 2024-11-25: Do not change the engine! RanecuEngine is the only one supported by our EventAction.
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  if (run_seed == -1)
    G4Random::setTheSeed(time(nullptr));
  else
    G4Random::setTheSeed(run_seed);

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive *ui = nullptr;
  if (!macro_path.size())
    ui = new G4UIExecutive(argc, argv, args["session"].as<G4String>());

  // Construct the default run manager
  //
  auto *runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  int threads = args["threads"].as<G4int>();
#ifdef G4MULTITHREADED
  if (threads > 0)
  {
    G4cout << " [GEANT4] setting: running with " << threads << " threads" << G4endl;
    runManager->SetNumberOfThreads(threads);
  }
#endif

  // Set mandatory initialization classes
  //
  auto detConstruction = new MuDetectorConstruction(args["detector"].as<G4String>(), args["export"].as<G4String>());
  runManager->SetUserInitialization(detConstruction);
  
  // Set physics list
  auto physicsList = new FTFP_BERT;
  // physicsList->RegisterPhysics(new ElectronEarthCutBuilder()); // Electron cuts
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  auto actionInitialization = new MuActionInitialization(detConstruction, args);
  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  //
  auto visManager = new G4VisExecutive("Quiet");
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Run in Batch mode / Interactive mode (Process macro or start UI session)
  //
  if (macro_path.size())
  {
    // batch mode
    // Setup aliases first (if there is any)
    for (std::size_t i = 1; i < macro_commands.size(); i += 2)
      UImanager->ApplyCommand("/control/alias " + macro_commands[i] + " " + macro_commands[i + 1]);
    UImanager->ApplyCommand("/control/execute " + macro_path);
  }
  else
  {
    // interactive mode : define UI session, with gui/vis setup macros
    UImanager->ApplyCommand("/control/execute " + PROJECT_SOURCE_DIR + "/macros/init_vis.mac");
    // if (ui->IsGUI())
    UImanager->ApplyCommand("/control/execute  " + PROJECT_SOURCE_DIR + "/macros/init_gui.mac");
    ui->SessionStart();
    delete ui;
  }

  // ** Interactive mode should never reach here **
  // Some special treatment in batch mode
  if (macro_path.size())
  {
    // For `recreate` and `filereader` generator in batch mode, run on all available events
    auto macro_string = util::io::readFileToString(macro_path);
    if ((util::notstd::strfind("recreate", macro_string) || util::notstd::strfind("filereader", macro_string)) && !util::notstd::strfind("/run/beamOn", macro_string))
    {

      runManager->Initialize();
      // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
      auto generator_name = MuPrimaryGeneratorAction::GetName();
      if (generator_name == "recreate" || generator_name == "filereader")
      {
        // const MuPrimaryGeneratorAction * generatoraction = static_cast<const MuPrimaryGeneratorAction*> (runManager->GetUserPrimaryGeneratorAction());
        // int numberOfEvent = generatoraction->GetGenerator()->GetEntries();
        int numberOfEvent = (MuPrimaryGeneratorAction::GetGenerator())->GetEntries();
        print(" [Generator] using recreate generator", generator_name, "with", numberOfEvent, "events");
        // runManager->BeamOn(numberOfEvent);
        UImanager->ApplyCommand("/run/beamOn " + std::to_string(numberOfEvent));
      }
    }
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
