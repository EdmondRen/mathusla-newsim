#ifndef MU__cry_hh
#define MU__cry_hh

// Geant4
#include "G4ParticleGun.hh"

// CRY
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"


// Project
#include "_Generator.hh"

namespace MuGenerators
{
    enum GEOTYPE
    {
        rect,
        box,
        sphere,
    };

    struct Vec3
    {
        double x, y, z;
    };        


    /// The primary generator action class with particle gum.
    ///
    /// It defines a single particle which hits the calorimeter
    /// perpendicular to the input face. The type of the particle
    /// can be changed via the G4 build-in commands of G4ParticleGun class
    /// (see the macros provided with this example).
    class MuCRY : public Generator
    {
    public:
        MuCRY(const std::string &name,
              const std::string &description,
              const std::string &PROJECT_SOURCE_DIR);
        virtual ~MuCRY() = default;

        // Core function 1: GeneratePrimaryVertex()
        // This will be called by GeneratorAction::GeneratePrimaries()
        void GeneratePrimaryVertex(G4Event *event) override;

        // Core function 2: GeneratePrimaryVertex()
        // This is used to set generator parameters
        void SetNewValue(G4UIcommand *command,
                         G4String value) override;
        float extractSubBoxLength(const std::string& filename);

        // Other helper functions
        // std::ostream &Print(std::ostream &os = std::cout) const override;

    private:



        G4ParticleGun *fParticleGun;
        G4ParticleTable *fparticleTable;
        
        CRYGenerator *fCRYgenerator;
        std::vector<CRYParticle *> *cry_generated;
        std::map<std::string, float> fCRY_additional_setup;
        
        std::string PROJECT_SOURCE_DIR;

        // Two corners of the box
        Vec3 subBoxMin,subBoxMax;
        GEOTYPE samplingShape;

        // Messenger commands
        G4UIcmdWithAString *_ui_pathname;
        G4UIcmdWith3VectorAndUnit *_ui_offset;
        G4UIcmdWithADoubleAndUnit *_ui_offset_t_low;
        G4UIcmdWithADoubleAndUnit *_ui_offset_t_high;
        G4UIcmdWithADoubleAndUnit *_ui_ekin_low;
        G4UIcmdWithADoubleAndUnit *_ui_ekin_high;
        G4UIcmdWithADouble *_ui_particle;
    };
}

#endif
