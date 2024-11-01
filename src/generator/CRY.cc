// GEANT4
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomTools.hh" // For G4UniformRand

// CRY
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"

// Project
#include "generator/CRY.hh"
#include "util.hh"

namespace MuGenerators
{

    //--------------------------------------------------------------------------
    // Wrapper for random number generator
    template <class T>
    class RNGWrapper
    {
    public:
        static void set(T *object, double (T::*func)(void));
        static double rng(void);

    private:
        static T *m_obj;
        static double (T::*m_func)(void);
    };

    template <class T>
    T *RNGWrapper<T>::m_obj;

    template <class T>
    double (T::*RNGWrapper<T>::m_func)(void);

    template <class T>
    void RNGWrapper<T>::set(T *object, double (T::*func)(void))
    {
        m_obj = object;
        m_func = func;
    }

    template <class T>
    double RNGWrapper<T>::rng(void) { return (m_obj->*m_func)(); }

    double GenerateRandomInRange(double min, double max)
    {
        if (min >= max)
        {
            throw std::invalid_argument("Invalid range: min must be less than max");
        }
        // Generate a random number in the range [min, max)
        return min + (max - min) * G4UniformRand();
    }
    //--------------------------------------------------------------------------

    MuCRY::MuCRY(const std::string &name,
                 const std::string &description,
                 const std::string &project_source_dir) : Generator(name, description), PROJECT_SOURCE_DIR(project_source_dir)
    {
        // A particle gun to use with CRY
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);
        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particleDefinition);
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 20 * m));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
        fParticleGun->SetParticleEnergy(50. * MeV);

        // CRY initialization
        auto cry_setupString = util::io::readFileToString_CRY(PROJECT_SOURCE_DIR + "/macros/generators/cry_default.file");

        G4cout << cry_setupString;

        CRYSetup *cry_setup = new CRYSetup(cry_setupString, PROJECT_SOURCE_DIR + "/cry_v1.7/data");
        // Set random number generator to use GEANT4 engine
        RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
        cry_setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
        // Make the CRY generator
        this->fCRYgenerator = new CRYGenerator(cry_setup);
        this->cry_generated = new std::vector<CRYParticle *>; // vector to hold generated particles

        // Create the table containing all particle names
        this->fparticleTable = G4ParticleTable::GetParticleTable();

        // Initialize Other parameters that is not part of CRY, but are needed to generate particle at correct points
        fCRY_additional_setup["offset_x"] = 0 * m;
        fCRY_additional_setup["offset_y"] = 0 * m;
        fCRY_additional_setup["offset_z"] = 4 * m;
        fCRY_additional_setup["offset_t_low"] = -1000 * ns;
        fCRY_additional_setup["offset_t_high"] = 1000 * ns;
        fCRY_additional_setup["ekin_cut_low"] = 0.1 * GeV;
        fCRY_additional_setup["ekin_cut_high"] = 100 * TeV;

        // Make messenger commands
        _ui_pathname = CreateCommand<G4UIcmdWithAString>("pathname", "Set pathname of CRY parameters file.");
        _ui_pathname->SetParameterName("pathname", false, false);
        _ui_pathname->AvailableForStates(G4State_PreInit, G4State_Idle);
        _ui_offset = CreateCommand<G4UIcmdWith3VectorAndUnit>("offset", "Set x-y-z offset with unit.");
        _ui_offset->SetParameterName("offset_x", "offset_y", "offset_z", false, false);
        _ui_offset->AvailableForStates(G4State_PreInit, G4State_Idle);    
        _ui_offset_t_low = CreateCommand<G4UIcmdWithADoubleAndUnit>("offset_t_low", "Set time offset range lower bound with unit.");
        _ui_offset_t_low->SetParameterName("offset_t_low", false, false);
        _ui_offset_t_low->AvailableForStates(G4State_PreInit, G4State_Idle);
        _ui_offset_t_high = CreateCommand<G4UIcmdWithADoubleAndUnit>("offset_t_high", "Set time offset range upper bound with unit.");
        _ui_offset_t_high->SetParameterName("offset_t_high", false, false);
        _ui_offset_t_high->AvailableForStates(G4State_PreInit, G4State_Idle);
        _ui_ekin_low = CreateCommand<G4UIcmdWithADoubleAndUnit>("ekin_low", "Set time kinetic energy cut lower bound with unit.");
        _ui_ekin_low->SetParameterName("ekin_low", false, false);
        _ui_ekin_low->AvailableForStates(G4State_PreInit, G4State_Idle);
        _ui_ekin_high = CreateCommand<G4UIcmdWithADoubleAndUnit>("ekin_high", "Set time kinetic energy cut upper bound with unit.");
        _ui_ekin_high->SetParameterName("ekin_high", false, false);
        _ui_ekin_high->AvailableForStates(G4State_PreInit, G4State_Idle);                                  
    }

    // Core function 1: GeneratePrimaryVertex()
    // This will be called by GeneratorAction::GeneratePrimaries()
    void MuCRY::GeneratePrimaryVertex(G4Event *anEvent)
    {

        G4String particleName;
        bool pass_cuts = false;

        int countAttempt = 0;
        do
        {
            cry_generated->clear();
            fCRYgenerator->genEvent(cry_generated);

            for (unsigned j = 0; j < cry_generated->size(); j++)
            {

                particleName = CRYUtils::partName((*cry_generated)[j]->id());
                G4ParticleDefinition *particleDefinition = fparticleTable->FindParticle((*cry_generated)[j]->PDGid());

                double kinEnergy = (*cry_generated)[j]->ke() * MeV;
                if (kinEnergy >= fCRY_additional_setup["ekin_cut_low"] && kinEnergy <= fCRY_additional_setup["ekin_cut_high"] )
                    pass_cuts = true; // 2670
                else
                    countAttempt++;
            }

        } while (pass_cuts == false);

        //....debug output
        G4cout << "\nEvent=" << anEvent->GetEventID() << " "
               << "CRY generated nparticles=" << cry_generated->size()
               << " pass Ekin threshold: " << pass_cuts << G4endl;

        // Sample a time for this event
        G4double t0 = GenerateRandomInRange(fCRY_additional_setup["offset_t_low"], fCRY_additional_setup["offset_t_high"]);

        for (unsigned j = 0; j < cry_generated->size(); j++)
        {
            particleName = CRYUtils::partName((*cry_generated)[j]->id());

            auto particleDefinition = fparticleTable->FindParticle((*cry_generated)[j]->PDGid());

            G4double fParticleEkin = (*cry_generated)[j]->ke() * MeV;
            G4double fParticleMass = particleDefinition->GetPDGMass() * MeV;
            G4double fParticleMomentum = sqrt(fParticleEkin * fParticleEkin + 2 * fParticleEkin * fParticleMass);
            G4double fParticlePosX = (*cry_generated)[j]->x() * m;
            G4double fParticlePosY = (*cry_generated)[j]->y() * m;
            G4double fParticlePosZ = (*cry_generated)[j]->z() * m;
            G4double fParticleMomentumDirectionU = (*cry_generated)[j]->u();
            G4double fParticleMomentumDirectionV = (*cry_generated)[j]->v();
            G4double fParticleMomentumDirectionW = (*cry_generated)[j]->w();
            G4double fParticleMomentumX = fParticleMomentum * fParticleMomentumDirectionU;
            G4double fParticleMomentumY = fParticleMomentum * fParticleMomentumDirectionV;
            G4double fParticleMomentumZ = fParticleMomentum * fParticleMomentumDirectionW;
            G4double fParticleTime = t0; //(*cry_generated)[j]->t();

            fParticleGun->SetParticleDefinition(particleDefinition);
            fParticleGun->SetParticleMomentum(sqrt(fParticleEkin * fParticleEkin + 2 * fParticleEkin * fParticleMass));
            fParticleGun->SetParticlePosition(G4ThreeVector(fParticlePosX + fCRY_additional_setup["offset_x"],
                                                            fParticlePosY + fCRY_additional_setup["offset_y"],
                                                            fParticlePosZ + fCRY_additional_setup["offset_z"]));
            fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fParticleMomentumDirectionU, fParticleMomentumDirectionV, fParticleMomentumDirectionW));
            fParticleGun->SetParticleTime(fParticleTime);
            fParticleGun->GeneratePrimaryVertex(anEvent);
        }
    }

    // Core function 2: GeneratePrimaryVertex()
    // This is used to set generator parameters
    void MuCRY::SetNewValue(G4UIcommand *command,
                            G4String value)
    {
        if (command == _ui_pathname)
        {
            // CRY initialization
            auto cry_setupString = util::io::readFileToString_CRY(value);
            G4cout << "\nCRY setup string: \n" <<cry_setupString << "\n";
            CRYSetup *cry_setup = new CRYSetup(cry_setupString, PROJECT_SOURCE_DIR + "/cry_v1.7/data");
            // Set random number generator to use GEANT4 engine
            RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
            cry_setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
            // Make the CRY generator
            this->fCRYgenerator = new CRYGenerator(cry_setup);
        }
        else if (command == _ui_offset)
        {
            auto vec = _ui_offset->GetNew3VectorValue(value);
            fCRY_additional_setup["offset_x"] = vec[0];
            fCRY_additional_setup["offset_y"] = vec[1];
            fCRY_additional_setup["offset_z"] = vec[2];            
        }
        else if (command == _ui_offset_t_low)
            fCRY_additional_setup["offset_t_low"] = _ui_offset_t_low->GetNewDoubleValue(value);
        else if (command == _ui_offset_t_high)
            fCRY_additional_setup["offset_t_high"] = _ui_offset_t_high->GetNewDoubleValue(value); 
        else if (command == _ui_ekin_low)
            fCRY_additional_setup["ekin_cut_low"] = _ui_ekin_low->GetNewDoubleValue(value);
        else if (command == _ui_ekin_high)
            fCRY_additional_setup["ekin_cut_high"] = _ui_ekin_high->GetNewDoubleValue(value);                        
    }

}