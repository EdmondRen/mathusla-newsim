#ifndef MU__recreate_hh
#define MU__recreate_hh

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
// Geant4
#include "G4Version.hh"
#if (G4VERSION_NUMBER < 1100)
#include "g4root.hh"
#else
#include "G4AnalysisManager.hh "
#endif

#include "G4RootAnalysisReader.hh"

// Project includes
#include "_Generator.hh"

namespace MuGenerators
{
    /// The primary generator action class for recreating previous events.
    class MuRecreate : public Generator
    {
    public:
        MuRecreate(const std::string &name,
                   const std::string &description,
                   const std::string &PROJECT_SOURCE_DIR);
        virtual ~MuRecreate() = default;

        // Core function 1: GeneratePrimaryVertex()
        // This will be called by GeneratorAction::GeneratePrimaries()
        void GeneratePrimaryVertex(G4Event *event) override;

        // Core function 2: GeneratePrimaryVertex()
        // This is used to set generator parameters
        void SetNewValue(G4UIcommand *command,
                         G4String value) override;

        // Other helper functions
        // std::ostream &Print(std::ostream &os = std::cout) const override;
        void ReadEventRecords(std::string filename);
        int GetEntries() const override;

    private:
        std::string PROJECT_SOURCE_DIR;

        int EVENTS_TOTAL, EVENTS_COUNTER;
        std::string root_filename;
        // G4AnalysisReader* analysisReader; // Tom: have to switch to pure ROOT approach.
        TTree *InputTree;
        TFile *InputFile;
        bool SetSeed;

        // Raw data holder
        int data_seed_init;
        int data_seed_0_raw;
        int data_seed_1_raw;
        std::vector<double> *data_pdgid = nullptr;
        std::vector<float> *data_x = nullptr;
        std::vector<float> *data_y = nullptr;
        std::vector<float> *data_z = nullptr;
        std::vector<float> *data_t = nullptr;
        std::vector<float> *data_px = nullptr;
        std::vector<float> *data_py = nullptr;
        std::vector<float> *data_pz = nullptr;
        std::vector<unsigned long> seed_combined;

        // Messenger commands
        G4UIcmdWithAString *_ui_pathname;
        G4UIcmdWithAnInteger *_ui_setentry;
    };
}

#endif