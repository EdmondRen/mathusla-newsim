#ifndef MU__GENERATOR_HH
#define MU__GENERATOR_HH
#pragma once

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

// G4 include
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <G4UIcommand.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UImessenger.hh>
#include <G4UImanager.hh>

namespace MuGenerators
{

    /*Default Particle Generator Template
    This is not directly related to Geant4.
    The basic requirement is a GeneratePrimaryVertex(G4Event*) function
    which add whatever particles to the G4Event
    */ 
    class Generator : public G4UImessenger
    {
    public:
        Generator(const std::string &name,
                  const std::string &description);

        virtual ~Generator() = default;

        // Core function 1: GeneratePrimaryVertex()
        // This will be called by GeneratorAction::GeneratePrimaries()
        virtual void GeneratePrimaryVertex(G4Event *event);

        // Core function 2: GeneratePrimaryVertex()
        // This is used to set generator parameters
        virtual void SetNewValue(G4UIcommand *command,
                                 G4String value) override;

        // Other helper functions
        virtual std::ostream &Print(std::ostream &os = std::cout) const;

        static const std::string MessengerDirectory;

    };


} // namespace Generator

#endif /* MU__PHYSICS_GENERATOR_HH */
