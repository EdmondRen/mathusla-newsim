

#include "generator/ParticleGun.hh"

namespace MuGenerators
{
    ParticleGun::ParticleGun(const std::string &name,
                             const std::string &description) : Generator(name, description)
    {
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        // Set the default parameters for particle gun
        // default particle kinematic
        //
        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particleDefinition);
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 20*m));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
        fParticleGun->SetParticleEnergy(50. * MeV);

        // In order to avoid dependence of PrimaryGeneratorAction
        // on DetectorConstruction class we get world volume
        // from G4LogicalVolumeStore
        //
        // G4double worldZHalfLength = 0.;
        // auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

        // // Check that the world volume has box shape
        // G4Box *worldBox = nullptr;
        // if (worldLV)
        //     worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());

        // if (worldBox)
        //     worldZHalfLength = worldBox->GetZHalfLength();
        // else
        // {
        //     G4ExceptionDescription msg;
        //     msg << "World volume of box shape not found." << G4endl;
        //     msg << "Perhaps you have changed geometry." << G4endl;
        //     msg << "The gun will be place in the center.";
        //     G4Exception("MuPrimaryGeneratorAction::GeneratePrimaries()",
        //                 "MyCode0002", JustWarning, msg);
        // }

        // // Set gun position to top at the middle of world. 
        // fParticleGun
        //     ->SetParticlePosition(G4ThreeVector(0., 0., worldZHalfLength));        

    }

    // Core function 1: GeneratePrimaryVertex()
    // This will be called by GeneratorAction::GeneratePrimaries()
    void ParticleGun::GeneratePrimaryVertex(G4Event *anEvent)
    {
        // Clear the store of generated particles
        this->genParticles.clear();


        // Run the generator
        fParticleGun->GeneratePrimaryVertex(anEvent);

        // Add the particle to store
        Particle newParticle = Particle(fParticleGun->GetParticleDefinition()->GetPDGEncoding(),
                                            fParticleGun->GetParticlePosition()[0],
                                            fParticleGun->GetParticlePosition()[1],
                                            fParticleGun->GetParticlePosition()[2],
                                            fParticleGun->GetParticleTime(),
                                            fParticleGun->GetParticleMomentum()*fParticleGun->GetParticleMomentumDirection()[0],
                                            fParticleGun->GetParticleMomentum()*fParticleGun->GetParticleMomentumDirection()[1],
                                            fParticleGun->GetParticleMomentum()*fParticleGun->GetParticleMomentumDirection()[2],
                                            0);
        this->genParticles.push_back(newParticle);        
    }

    // Core function 2: GeneratePrimaryVertex()
    // This is used to set generator parameters
    void ParticleGun::SetNewValue(G4UIcommand *command,
                                  G4String value)
    {
        (void)command;
        (void)value;
    }

    // Other helper functions
    std::ostream &ParticleGun::Print(std::ostream &os) const
    {
        return os;
    }

}