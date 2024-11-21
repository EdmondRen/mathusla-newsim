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
    // ----------------------------------------------------------------------------
    // Defining datatypes that are useful for generators
    // ----------------------------------------------------------------------------

    // Pseudo-Lorentz Invariant Triplet struct {pT, eta, phi}
    struct PseudoLorentzTriplet
    {
        double pT, eta, phi;
    };
    constexpr bool operator==(const PseudoLorentzTriplet &left,
                              const PseudoLorentzTriplet &right)
    {
        return left.pT == right.pT && left.eta == right.eta && left.phi == right.phi;
    }

    // Particle struct
    struct Particle
    {
        // Necessary parameters
        int pdgid; // PDG identifier
        double px, py, pz;
        double x, y, z, t;

        // Optional parameters
        int index; // A counter, can be initialized to the value you want

        Particle() = default;

        Particle(int identifier,
                 double x_momentum,
                 double y_momentum,
                 double z_momentum)
            : pdgid(identifier), px(x_momentum), py(y_momentum), pz(z_momentum), x(0), y(0), z(0), t(0) {}

        Particle(int identifier,
                 double x_position,
                 double y_position,
                 double z_position,
                 double t_position,
                 double x_momentum,
                 double y_momentum,
                 double z_momentum,
                 int genParticleIndex)
            : pdgid(identifier), px(x_momentum), py(y_momentum), pz(z_momentum),
              x(x_position), y(y_position), z(z_position), t(t_position), index(genParticleIndex) {}
    };

    // List of Particles
    using ParticleVector = std::vector<Particle>;

    // Vertex:
    struct Vertex
    {
        double tv, xv, yv, zv;
        ParticleVector particles;
    };

    // ----------------------------------------------------------------------------

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

        // Data Store: save the list of generated particles in the current event
        ParticleVector genParticles;

        // Other helper functions
        virtual std::ostream &Print(std::ostream &os = std::cout) const;

        static const std::string MessengerDirectory;
    };



} // namespace Generator

#endif /* MU__PHYSICS_GENERATOR_HH */
