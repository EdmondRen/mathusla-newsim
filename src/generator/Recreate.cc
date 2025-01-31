
// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
// Geant4
#include "Randomize.hh"
// Project
#include "generator/Recreate.hh"
#include "util.hh"

namespace MuGenerators
{
    MuRecreate::MuRecreate(const std::string &name,
                           const std::string &description,
                           const std::string &project_source_dir) : Generator(name, description), PROJECT_SOURCE_DIR(project_source_dir)
    {
        // Create variables
        // this->analysisReader = G4AnalysisReader::Instance();
        this->EVENTS_TOTAL = 0;
        this->EVENTS_COUNTER = 0;

        // Make messenger commands
        this->_ui_pathname = CreateCommand<G4UIcmdWithAString>("pathname", "Set pathname of PARMA configuration file.");
        this->_ui_pathname->SetParameterName("pathname", false, false);
        this->_ui_pathname->AvailableForStates(G4State_PreInit, G4State_Idle);
    }

    void MuRecreate::GeneratePrimaryVertex(G4Event *event)
    {
        // Clear the store of generated particles
        this->genParticles.clear();
        this->seed_combined.clear();
        // print("*************Total Entries", EVENTS_TOTAL);

        // if (this->analysisReader->GetNtupleRow())
        this->InputTree->GetEntry(this->EVENTS_COUNTER);
        if (EVENTS_COUNTER < EVENTS_TOTAL)
        {
            // Generate all particles in this event
            for (u_int i = 0; i < (*data_x).size(); i++)
            {
                // make a particle, then add it to the event and the particle store.
                Particle newParticle = Particle(static_cast<u_int64_t>((*data_pdgid)[i]),
                                                (*data_x)[i],
                                                (*data_y)[i],
                                                (*data_z)[i],
                                                (*data_t)[i],
                                                (*data_px)[i],
                                                (*data_py)[i],
                                                (*data_pz)[i],
                                                0);
                // Generate the particle
                const auto vertex = new G4PrimaryVertex(newParticle.x, newParticle.y, newParticle.z, newParticle.t);
                vertex->SetPrimary(new G4PrimaryParticle(newParticle.pdgid, newParticle.px, newParticle.py, newParticle.pz));
                event->AddPrimaryVertex(vertex);

                // Save all generated particles of the current event
                this->genParticles.push_back(newParticle);
            }

            if (SetSeed)
            {
                // Restore the random number generator status
                //  * make a state vector for random engine. Need to cast from int to unsigned int
                //  * RanecuEngine state has 4 values {pointer, theSeed, seed0, seed1}
                //  Two steps: set the seed (to update the internal seed table), the set the state
                seed_combined.push_back(0);
                seed_combined.push_back(1);
                seed_combined.push_back(data_seed_0_raw);
                seed_combined.push_back(data_seed_1_raw);
                CLHEP::HepRandom::getTheEngine()->getState(seed_combined);
            }
        }
        else
        {
            G4cout << " ERROR [Generator: Recreate] Reached the end of the input file, no more events to generate. Please change the number in /run/beamOn X to match the number of entries of the input ROOT file." << G4endl;
            exit(0);
        }

        this->EVENTS_COUNTER += 1;
    }

    int MuRecreate::GetEntries() const
    {
        return EVENTS_TOTAL;
    }

    void MuRecreate::SetNewValue(G4UIcommand *command,
                                 G4String value)
    {
        if (command == _ui_pathname)
        {
            this->root_filename = value;

            // Open the root file with analysis reader
            // Tom: won't work with Geant4.10. GetNtuple is badly implemented
            // analysisReader->SetFileName(value);
            // G4int ntupleid = analysisReader->GetNtuple("data");
            // auto tree = analysisReader->GetNtuple(ntupleid);
            // EVENTS_TOTAL = tree->entries();
            // analysisReader->SetNtupleIColumn("Seed_0", data_seed_0_raw);
            // analysisReader->SetNtupleIColumn("Seed_1", data_seed_1_raw);
            // analysisReader->SetNtupleDColumn("Gen_pdgID", data_pdgid);
            // analysisReader->SetNtupleFColumn("Gen_x", data_x);
            // analysisReader->SetNtupleFColumn("Gen_y", data_y);
            // analysisReader->SetNtupleFColumn("Gen_z", data_z);
            // analysisReader->SetNtupleFColumn("Gen_t", data_t);
            // analysisReader->SetNtupleFColumn("Gen_px", data_px);
            // analysisReader->SetNtupleFColumn("Gen_py", data_py);
            // analysisReader->SetNtupleFColumn("Gen_pz", data_pz);

            // Pure ROOT approach
            InputFile = TFile::Open(this->root_filename.c_str());
            if (!InputFile)
                return;

            auto input_tree_name = "data";
            InputTree = (TTree *)InputFile->Get(input_tree_name);
            if (!InputTree)
                return;

            this->EVENTS_TOTAL = InputTree->GetEntries();

            // Check if the Seed_0 branch exist:
            auto seed_found = InputTree->GetListOfBranches()->FindObject("Seed_0");
            if (seed_found == 0)
            {
                print("No seed info. Will create events without modifying random seeds");
                SetSeed = false;
            }
            else
            {
                SetSeed = true;
                InputTree->SetBranchAddress("Seed_0", &data_seed_0_raw);
                InputTree->SetBranchAddress("Seed_1", &data_seed_1_raw);
            }

            InputTree->SetBranchAddress("Gen_pdgID", &data_pdgid);
            InputTree->SetBranchAddress("Gen_x", &data_x);
            InputTree->SetBranchAddress("Gen_y", &data_y);
            InputTree->SetBranchAddress("Gen_z", &data_z);
            InputTree->SetBranchAddress("Gen_t", &data_t);
            InputTree->SetBranchAddress("Gen_px", &data_px);
            InputTree->SetBranchAddress("Gen_py", &data_py);
            InputTree->SetBranchAddress("Gen_pz", &data_pz);
        }
    }

}
