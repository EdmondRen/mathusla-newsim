
#include "Randomize.hh"

#include "Recreate.hh"

namespace MuGenerators
{
    MuRecreate::MuRecreate(const std::string &name,
                           const std::string &description,
                           const std::string &project_source_dir) : Generator(name, description), PROJECT_SOURCE_DIR(project_source_dir)
    {
        // Create variables
        this->analysisReader = G4AnalysisReader::Instance();
        this->EVENTS_TOTAL = 0;

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
        this->EVENTS_COUNTER += 1;

        if (this->analysisReader->GetNtupleRow())
        {
            for (int i=0; i<data_x.size(); i++)
            {
                // make a particle, then add it to the event and the particle store.
                Particle newParticle = Particle(static_cast<u_int64_t>(data_pdgid[i]),
                                                data_x[i],
                                                data_y[i],
                                                data_z[i],
                                                data_t[i],
                                                data_px[i],
                                                data_py[i],
                                                data_pz[i],
                                                0);

                // make a seed pair. Need to cast from int to unsigned int
                seed_combined.push_back(*reinterpret_cast<unsigned long *>(&data_seed_0_raw));
                seed_combined.push_back(*reinterpret_cast<unsigned long *>(&data_seed_1_raw));


                // Generate the particle
                const auto vertex = new G4PrimaryVertex(newParticle.x, newParticle.y, newParticle.z, newParticle.t);
                vertex->SetPrimary(new G4PrimaryParticle(newParticle.pdgid, newParticle.px, newParticle.py, newParticle.pz));
                event->AddPrimaryVertex(vertex);

                // Save all generated particles of the current event
                this->genParticles.push_back(newParticle);

                // Restore the random number generator status
                CLHEP::HepRandom::getTheEngine()->get(seed_combined);
            }
        }
        else
            G4cout<< " ERROR [Generator: Recreate] Reached the end of the input file, no more events to generate. Please change the number in /run/beamOn X to match the number of entries of the input ROOT file.";
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
            analysisReader->SetFileName(value);
            G4int ntupleid = analysisReader->GetNtuple("raw");
            auto tree = analysisReader->GetNtuple(ntupleid1234); 134
            EVENTS_TOTAL = tree->entries();
            analysisReader->SetNtupleIColumn("Seed_0", data_seed_0_raw);
            analysisReader->SetNtupleIColumn("Seed_1", data_seed_1_raw);
            analysisReader->SetNtupleDColumn("Gen_pdgID", data_pdgid);            
            analysisReader->SetNtupleFColumn("Gen_x", data_x);
            analysisReader->SetNtupleFColumn("Gen_y", data_y);
            analysisReader->SetNtupleFColumn("Gen_z", data_z);
            analysisReader->SetNtupleFColumn("Gen_t", data_t);
            analysisReader->SetNtupleFColumn("Gen_px", data_px);
            analysisReader->SetNtupleFColumn("Gen_py", data_py);
            analysisReader->SetNtupleFColumn("Gen_pz", data_pz);
        }
    }

}
