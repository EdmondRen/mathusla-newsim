#ifndef MU__parma_hh
#define MU__parma_hh

// Geant4
#include "G4ParticleGun.hh"

// Project includes
#include "_Generator.hh"
#include "generator/parma/parma_util.hh"

namespace MuGenerators
{
    /// The primary generator action class with PARMA.
    class MuPARMA : public Generator
    {
    public:
        MuPARMA(const std::string &name,
              const std::string &description,
              const std::string &PROJECT_SOURCE_DIR);

        // Core function 1: GeneratePrimaryVertex()
        // This will be called by GeneratorAction::GeneratePrimaries()
        void GeneratePrimaryVertex(G4Event *event) override;

        // Core function 2: GeneratePrimaryVertex()
        // This is used to set generator parameters
        void SetNewValue(G4UIcommand *command,
                         G4String value) override;

        // Get methods
        std::map<std::string, std::string> getMetaData() override;
        float getEventWeight() override;

        // Other helper functions
        // std::ostream &Print(std::ostream &os = std::cout) const override;
        void startPARMA(const std::string &config_filename);
        void resetDimensions();

    private:

        G4ParticleTable *fparticleTable;
        G4ParticleGun *fparticleGun;
        
        PARMA::ParmaGen fPARMAgenerator;
        PARMA::ParmaParticle parma_generated;
        std::map<std::string, float> fPARMA_additional_config;
        
        std::string PROJECT_SOURCE_DIR;

        // Event information
        long event_counter;
        float event_weight;
        std::map<std::string, std::string> metadata;

        // Two corners of the box
        double subboxLength;
        Vec3 subBoxMin,subBoxMax;
        GEOTYPE samplingShape;

        // Messenger commands
        G4UIcmdWithAString *_ui_pathname;
        G4UIcmdWithAnInteger *_ui_shape;
        G4UIcmdWith3VectorAndUnit *_ui_box;
        G4UIcmdWith3VectorAndUnit *_ui_offset;
        G4UIcmdWithADoubleAndUnit *_ui_offset_t_low;
        G4UIcmdWithADoubleAndUnit *_ui_offset_t_high;
        G4UIcmdWithADoubleAndUnit *_ui_ekin_low;
        G4UIcmdWithADoubleAndUnit *_ui_ekin_high;
        G4UIcmdWithADouble *_ui_particle;
        G4UIcommand *_ui_update;
    };
}

#endif
