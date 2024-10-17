
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "Generator.hh"

namespace MuGenerators
{   
    /// The primary generator action class with particle gum.
    ///
    /// It defines a single particle which hits the calorimeter 
    /// perpendicular to the input face. The type of the particle
    /// can be changed via the G4 build-in commands of G4ParticleGun class 
    /// (see the macros provided with this example).
    class ParticleGun : public Generator
    {
    public:
        ParticleGun(const std::string &name,
                    const std::string &description);
        virtual ~ParticleGun() = default;

        // Core function 1: GeneratePrimaryVertex()
        // This will be called by GeneratorAction::GeneratePrimaries()
        void GeneratePrimaryVertex(G4Event *event) override;

        // Core function 2: GeneratePrimaryVertex()
        // This is used to set generator parameters
        void SetNewValue(G4UIcommand *command,
                         G4String value) override;

        // Other helper functions
        std::ostream &Print(std::ostream &os = std::cout) const override;
    private:
        G4ParticleGun* fParticleGun;
    };
}
