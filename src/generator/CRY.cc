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

    // Helper function to find intersection range
    bool intersectSlab(double p0, double d, double min, double max, double &tmin, double &tmax)
    {
        if (std::abs(d) < 1e-8)
        {
            // Line is parallel to the slab
            return p0 >= min && p0 <= max;
        }
        double t1 = (min - p0) / d;
        double t2 = (max - p0) / d;
        if (t1 > t2)
            std::swap(t1, t2);
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        return tmin <= tmax;
    }

    bool doesLineIntersectBox(const Vec3 &p0, const Vec3 &d, const Vec3 &boxMin, const Vec3 &boxMax)
    {
        double tmin = -INFINITY, tmax = INFINITY;

        // Check x-axis slab
        if (!intersectSlab(p0.x, d.x, boxMin.x, boxMax.x, tmin, tmax))
            return false;

        // Check y-axis slab
        if (!intersectSlab(p0.y, d.y, boxMin.y, boxMax.y, tmin, tmax))
            return false;

        // Check z-axis slab
        if (!intersectSlab(p0.z, d.z, boxMin.z, boxMax.z, tmin, tmax))
            return false;

        return true;
    }

    MuCRY::MuCRY(const std::string &name,
                 const std::string &description,
                 const std::string &project_source_dir) : Generator(name, description), PROJECT_SOURCE_DIR(project_source_dir)
    {
        // Initialize Other parameters that is not part of CRY, but are needed to generate particle at correct points
        fCRY_additional_setup["use_shape"] = 0; // 0: rect, 1: box, 2: sphere

        fCRY_additional_setup["rect_lenx"] = 40 * m;
        fCRY_additional_setup["rect_leny"] = 4 * m;
        fCRY_additional_setup["rect_lenz"] = 10 * m;
        fCRY_additional_setup["sphere_r"] = 1 * m;

        fCRY_additional_setup["offset_x"] = 0 * m;
        fCRY_additional_setup["offset_y"] = 0 * m;
        fCRY_additional_setup["offset_z"] = fCRY_additional_setup["rect_lenz"];
        fCRY_additional_setup["offset_t_low"] = -1000 * ns;
        fCRY_additional_setup["offset_t_high"] = 1000 * ns;

        fCRY_additional_setup["ekin_cut_low"] = 0.1 * GeV;
        fCRY_additional_setup["ekin_cut_high"] = 100 * TeV;
        fCRY_additional_setup["particle_pdgid"] = 0;

        // A particle gun to use with CRY
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);
        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particleDefinition);
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 20 * m));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
        fParticleGun->SetParticleEnergy(50. * MeV);

        // CRY initialization, using default cry config file
        auto cry_defaultConfig = "/macros/generators/cry_default.file";
        auto cry_setupString = util::io::readFileToString_CRY(PROJECT_SOURCE_DIR + cry_defaultConfig);
        // G4cout << cry_setupString;
        CRYSetup *cry_setup = new CRYSetup(cry_setupString, PROJECT_SOURCE_DIR + "/third_party/cry_v1.7/data");

        // Find the box width / sphere radius from the cry config file
        this->samplingShape = static_cast<GEOTYPE>(static_cast<int>(fCRY_additional_setup["use_shape"]));
        this->subBoxMin = {fCRY_additional_setup["offset_x"] - fCRY_additional_setup["rect_lenx"] / 2,
                           fCRY_additional_setup["offset_y"] - fCRY_additional_setup["rect_leny"] / 2,
                           fCRY_additional_setup["offset_z"] - fCRY_additional_setup["rect_lenz"]};
        this->subBoxMax = {fCRY_additional_setup["offset_x"] + fCRY_additional_setup["rect_lenx"] / 2,
                           fCRY_additional_setup["offset_y"] + fCRY_additional_setup["rect_leny"] / 2,
                           fCRY_additional_setup["offset_z"]};

        // Set random number generator to use GEANT4 engine
        RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
        cry_setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);

        // Make the CRY generator
        this->fCRYgenerator = new CRYGenerator(cry_setup);
        this->cry_generated = new std::vector<CRYParticle *>; // vector to hold generated particles

        // Create the table containing all particle names
        this->fparticleTable = G4ParticleTable::GetParticleTable();

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
        _ui_particle = CreateCommand<G4UIcmdWithADouble>("particle_pdgid", "Select the PDG ID of the particle to apply the cut on. Default is all particles. If this command is set, then only events CONTAINING the selected particle will be kept. To undo it, set the value to 0");
        _ui_particle->SetParameterName("particle_pdgid", false, false);
        _ui_particle->AvailableForStates(G4State_PreInit, G4State_Idle);
    }

    // Core function 1: GeneratePrimaryVertex()
    // This will be called by GeneratorAction::GeneratePrimaries()
    void MuCRY::GeneratePrimaryVertex(G4Event *anEvent)
    {
        // Clear the store of generated particles
        this->genParticles.clear();

        G4String particleName;
        double kinEnergy = 0;
        int pdgid = 0;
        bool pass_cuts1 = false;
        bool pass_cuts2 = false;

        int countAttempt = 0;
        double tmin = 1e100;
        do
        {
            cry_generated->clear();
            fCRYgenerator->genEvent(cry_generated);

            for (unsigned j = 0; j < cry_generated->size(); j++)
            {
                // G4ParticleDefinition *particleDefinition = fparticleTable->FindParticle((*cry_generated)[j]->PDGid());
                // particleName = CRYUtils::partName((*cry_generated)[j]->id());
                kinEnergy = (*cry_generated)[j]->ke() * MeV;
                pdgid = (*cry_generated)[j]->PDGid();

                // Cuts
                // 1. Select particles in the given energy range and with certain PDG id
                // If "particle_pdgid" is set, cut on that as well.
                if ((fCRY_additional_setup["particle_pdgid"] == 0 || fCRY_additional_setup["particle_pdgid"] == pdgid) && (kinEnergy >= fCRY_additional_setup["ekin_cut_low"] && kinEnergy <= fCRY_additional_setup["ekin_cut_high"]))
                    pass_cuts1 = true;
                else
                    countAttempt++;

                // 2. Only select events with tracks that intersect with the box
                if (this->samplingShape == box)
                {
                    G4double fParticlePosX = (*cry_generated)[j]->x() * m + fCRY_additional_setup["offset_x"];
                    G4double fParticlePosY = (*cry_generated)[j]->y() * m + fCRY_additional_setup["offset_y"];
                    G4double fParticlePosZ = (*cry_generated)[j]->z() * m + fCRY_additional_setup["offset_z"];
                    G4double fParticlePu = (*cry_generated)[j]->u();
                    G4double fParticlePv = (*cry_generated)[j]->v();
                    G4double fParticlePw = (*cry_generated)[j]->w();

                    Vec3 line_point = {fParticlePosX, fParticlePosY, fParticlePosZ};
                    Vec3 line_direction = {fParticlePu, fParticlePv, fParticlePw};
                    if (((*cry_generated)[j]->x() * m < ) || doesLineIntersectBox(line_point, line_direction, this->subBoxMin, this->subBoxMax))
                        pass_cuts2 = true;
                }
                else
                    pass_cuts2=true;

                // Find the time of the first particle
                if ((*cry_generated)[j]->t() < tmin)
                    tmin = (*cry_generated)[j]->t();
            }

        } while (pass_cuts1 == false || pass_cuts2 == false );

        //....debug output
        // G4cout << "\nEvent=" << anEvent->GetEventID() << " "
        //        << "CRY generated nparticles=" << cry_generated->size()
        //        << " pass Ekin threshold: " << pass_cuts << G4endl;
        // util::py::print(cry_generated->size());

        // Sample a time for this event
        G4double t0 = GenerateRandomInRange(fCRY_additional_setup["offset_t_low"], fCRY_additional_setup["offset_t_high"]);

        for (unsigned j = 0; j < cry_generated->size(); j++)
        {
            particleName = CRYUtils::partName((*cry_generated)[j]->id());

            auto pdgID = (*cry_generated)[j]->PDGid();
            auto particleDefinition = fparticleTable->FindParticle(pdgID);

            G4double fParticleEkin = (*cry_generated)[j]->ke() * MeV;
            G4double fParticleMass = particleDefinition->GetPDGMass() * MeV;
            G4double fParticlePosX = (*cry_generated)[j]->x() * m + fCRY_additional_setup["offset_x"];
            G4double fParticlePosY = (*cry_generated)[j]->y() * m + fCRY_additional_setup["offset_y"];
            G4double fParticlePosZ = (*cry_generated)[j]->z() * m + fCRY_additional_setup["offset_z"];
            G4double fParticleMomentumDirectionU = (*cry_generated)[j]->u();
            G4double fParticleMomentumDirectionV = (*cry_generated)[j]->v();
            G4double fParticleMomentumDirectionW = (*cry_generated)[j]->w();
            G4double fParticleMomentum = sqrt(fParticleEkin * fParticleEkin + 2 * fParticleEkin * fParticleMass);
            G4double fParticleMomentumX = fParticleMomentum * fParticleMomentumDirectionU;
            G4double fParticleMomentumY = fParticleMomentum * fParticleMomentumDirectionV;
            G4double fParticleMomentumZ = fParticleMomentum * fParticleMomentumDirectionW;
            G4double fParticleTime = t0 + ((*cry_generated)[j]->t() - tmin);

            // Generate the particle.
            // You can continue using particle gun,
            // fParticleGun->SetParticleDefinition(particleDefinition);
            // fParticleGun->SetParticlePosition(G4ThreeVector(fParticlePosX, fParticlePosY, fParticlePosZ));
            // fParticleGun->SetParticleMomentum(G4ThreeVector(fParticleMomentumX, fParticleMomentumY, fParticleMomentumZ));
            // fParticleGun->SetParticleTime(fParticleTime);
            // fParticleGun->GeneratePrimaryVertex(anEvent);

            // or you can make a particle, then add it to the event and the particle store.
            Particle newParticle = Particle(pdgID,
                                            fParticlePosX,
                                            fParticlePosY,
                                            fParticlePosZ,
                                            fParticleTime,
                                            fParticleMomentumX,
                                            fParticleMomentumY,
                                            fParticleMomentumZ,
                                            j);

            const auto vertex = new G4PrimaryVertex(newParticle.x, newParticle.y, newParticle.z, newParticle.t);
            vertex->SetPrimary(new G4PrimaryParticle(newParticle.pdgid, newParticle.px, newParticle.py, newParticle.pz));
            anEvent->AddPrimaryVertex(vertex);

            genParticles.push_back(newParticle);
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
            G4cout << "\nCRY setup string: \n"
                   << cry_setupString << "\n";
            CRYSetup *cry_setup = new CRYSetup(cry_setupString, PROJECT_SOURCE_DIR + "/third_party/cry_v1.7/data");
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
        else if (command == _ui_particle)
            fCRY_additional_setup["particle_pdgid"] = _ui_particle->GetNewDoubleValue(value);
    }

    float MuCRY::extractSubBoxLength(const std::string &filename)
    {
        std::ifstream file(filename); // Open the file
        if (!file.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return -1; // Return -1 to indicate an error
        }

        std::string line;
        while (std::getline(file, line))
        { // Read the file line by line
            std::istringstream iss(line);
            std::string key;
            int value;
            iss >> key >> value; // Extract key and value
            if (key == "subboxLength")
            {                 // Check if the key matches
                return value; // Return the extracted value
            }
        }

        std::cerr << "Error: 'subboxLength' not found in the file." << std::endl;
        return -1; // Return -1 to indicate the key was not found
    }
}