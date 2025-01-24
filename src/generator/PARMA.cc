#include <cmath>

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

// Project
#include "generator/PARMA.hh"
#include "generator/parma/parma_util.hh"
#include "util.hh"

namespace MuGenerators
{
    MuPARMA::MuPARMA(const std::string &name,
                     const std::string &description,
                     const std::string &project_source_dir) : Generator(name, description), PROJECT_SOURCE_DIR(project_source_dir)
    {
        // Initialize Other parameters that is not part of PARMA, but are needed to generate particle at correct points
        fPARMA_additional_config["use_shape"] = 0; // 0: square, 1: box, 2: sphere

        fPARMA_additional_config["box_lenx"] = 2 * m;
        fPARMA_additional_config["box_leny"] = 2 * m;
        fPARMA_additional_config["box_lenz"] = 1 * m;
        fPARMA_additional_config["sphere_r"] = 1 * m;

        fPARMA_additional_config["offset_x"] = 0 * m;
        fPARMA_additional_config["offset_y"] = 0 * m;
        fPARMA_additional_config["offset_z"] = fPARMA_additional_config["box_lenz"];
        fPARMA_additional_config["offset_t_low"] = -1000 * ns;
        fPARMA_additional_config["offset_t_high"] = 1000 * ns;

        fPARMA_additional_config["ekin_cut_low"] = 0.1 * GeV;
        fPARMA_additional_config["ekin_cut_high"] = 100 * TeV;
        fPARMA_additional_config["particle_pdgid"] = 0;

        fPARMA_additional_config["weight_function"] = 0; // 0: none, 1: weight1(), 2: weight2()

        // PARMA initialization, using default  config file
        auto parma_defaultConfig = "/macros/generators/parma_default.conf";
        G4cout<< "Start PARMA using configuration file: "<< PROJECT_SOURCE_DIR + parma_defaultConfig << G4endl;

        // Instantiate and config the PARMA generator
        this->fPARMAgenerator = PARMA::ParmaGen();

        // Set random number generator to use GEANT4 engine
        RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
        fPARMAgenerator.setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);

        startPARMA(PROJECT_SOURCE_DIR + parma_defaultConfig);

        // Create the table containing all particle names
        this->fparticleTable = G4ParticleTable::GetParticleTable();

        // Make messenger commands
        _ui_pathname = CreateCommand<G4UIcmdWithAString>("pathname", "Set pathname of PARMA configuration file.");
        _ui_pathname->SetParameterName("pathname", false, false);
        _ui_pathname->AvailableForStates(G4State_PreInit, G4State_Idle);
        _ui_shape = CreateCommand<G4UIcmdWithAnInteger>("shape", "Set the shape of generator. 0: square, 1: box, 2: sphere (dome)");
        _ui_shape->SetParameterName("shape", false, false);
        _ui_shape->AvailableForStates(G4State_PreInit, G4State_Idle);
        _ui_box = CreateCommand<G4UIcmdWith3VectorAndUnit>("box", "Set x-y-z length of the box.");
        _ui_box->SetParameterName("box_lenx", "box_leny", "box_lenz", false, false);
        _ui_box->AvailableForStates(G4State_PreInit, G4State_Idle);
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
        _ui_update = CreateCommand<G4UIcommand>("update", "Update the PARMA generator settings.");
        _ui_update->AvailableForStates(G4State_PreInit, G4State_Idle);        
    }

    void MuPARMA::startPARMA(const std::string &config_filename)
    {
        // Get the configuration file
        auto parcard = util::io::ParHandler(config_filename);
        auto config = parcard.GetConfig();


        // Setup the generator with the config map
        fPARMAgenerator.configure(config);
        subboxLength = fPARMAgenerator.subboxlength * cm;

        // Debug print
        // print("Parma config 0",
        //         fPARMAgenerator.amin,
        //         fPARMAgenerator.amax,
        //         fPARMAgenerator.emin,
        //         fPARMAgenerator.emax,
        //         fPARMAgenerator.alti,
        //         fPARMAgenerator.glat,
        //         fPARMAgenerator.glong,
        //         fPARMAgenerator.ip,
        //         fPARMAgenerator.iyear,
        //         fPARMAgenerator.imonth,
        //         fPARMAgenerator.iday
        //         );
    }

    void MuPARMA::resetDimensions()
    {
        this->samplingShape = static_cast<GEOTYPE>(static_cast<int>(fPARMA_additional_config["use_shape"]));
        if (samplingShape == box &&
            (subboxLength < fPARMA_additional_config["box_lenx"] ||
             subboxLength < fPARMA_additional_config["box_leny"] ))
        {
            print(subboxLength, fPARMA_additional_config["box_lenx"]);
            print("ERROR [Generator PARMA]: subboxLength in the config file does not cover the requested box length in GEANT4");
            exit(0);
        }
        this->subBoxMin = {fPARMA_additional_config["offset_x"] - fPARMA_additional_config["box_lenx"] / 2,
                           fPARMA_additional_config["offset_y"] - fPARMA_additional_config["box_leny"] / 2,
                           fPARMA_additional_config["offset_z"] - fPARMA_additional_config["box_lenz"]};
        this->subBoxMax = {fPARMA_additional_config["offset_x"] + fPARMA_additional_config["box_lenx"] / 2,
                           fPARMA_additional_config["offset_y"] + fPARMA_additional_config["box_leny"] / 2,
                           fPARMA_additional_config["offset_z"]};
    }

    // Core function 1: GeneratePrimaryVertex()
    // This will be called by GeneratorAction::GeneratePrimaries()
    void MuPARMA::GeneratePrimaryVertex(G4Event *anEvent)
    {
        // Clear the store of generated particles
        this->genParticles.clear();

        this->event_counter +=1;

        // Generate one event
        parma_generated = fPARMAgenerator.Generate();

        // Reject events whose tracks do not intersect with the box
        bool pass_cuts2 = false;
        if (this->samplingShape == box)
        {
            G4double fParticlePosX = parma_generated.x * cm + fPARMA_additional_config["offset_x"];
            G4double fParticlePosY = parma_generated.y * cm + fPARMA_additional_config["offset_y"];
            G4double fParticlePosZ = parma_generated.z * cm + fPARMA_additional_config["offset_z"];
            G4double fParticlePu = parma_generated.u;
            G4double fParticlePv = parma_generated.v;
            G4double fParticlePw = parma_generated.w;

            Vec3 line_point = {fParticlePosX, fParticlePosY, fParticlePosZ};
            Vec3 line_direction = {fParticlePu, fParticlePv, fParticlePw};
            if ((std::abs(parma_generated.x * cm) < fPARMA_additional_config["box_lenx"] * 0.5 &&
                 std::abs(parma_generated.y * cm) < fPARMA_additional_config["box_leny"] * 0.5) ||
                doesLineIntersectBox(line_point, line_direction, this->subBoxMin, this->subBoxMax))
            {
                pass_cuts2 = true;
            }
            else{
                // // Some debug information
                // print(subBoxMin.x,subBoxMin.y,subBoxMin.z);
                // print(subBoxMax.x,subBoxMax.y,subBoxMax.z);
                // print(line_point.x,line_point.y,line_point.z);
                // print(line_direction.x,line_direction.y,line_direction.z);
                // print("failed");
            }
        }
        else
            pass_cuts2 = true;
        if (!pass_cuts2)
            return;


        // Sample a time for this event
        G4double t0 = GenerateRandomInRange(fPARMA_additional_config["offset_t_low"], fPARMA_additional_config["offset_t_high"]);

        auto pdgID = parma_generated.pdgid;


        auto particleDefinition = fparticleTable->FindParticle(pdgID);

        G4double fParticleEkin = parma_generated.ke * MeV;
        G4double fParticleMass = particleDefinition->GetPDGMass() * MeV;
        G4double fParticlePosX = parma_generated.x * cm + fPARMA_additional_config["offset_x"];
        G4double fParticlePosY = parma_generated.y * cm + fPARMA_additional_config["offset_y"];
        G4double fParticlePosZ = parma_generated.z * cm + fPARMA_additional_config["offset_z"];
        G4double fParticleMomentumDirectionU = parma_generated.u;
        G4double fParticleMomentumDirectionV = parma_generated.v;
        G4double fParticleMomentumDirectionW = parma_generated.w;
        G4double fParticleMomentum = sqrt(fParticleEkin * fParticleEkin + 2 * fParticleEkin * fParticleMass);
        G4double fParticleMomentumX = fParticleMomentum * fParticleMomentumDirectionU;
        G4double fParticleMomentumY = fParticleMomentum * fParticleMomentumDirectionV;
        G4double fParticleMomentumZ = fParticleMomentum * fParticleMomentumDirectionW;
        G4double fParticleTime = t0 + parma_generated.t;

        // or you can make a particle, then add it to the event and the particle store.
        Particle newParticle = Particle(pdgID,
                                        fParticlePosX,
                                        fParticlePosY,
                                        fParticlePosZ,
                                        fParticleTime,
                                        fParticleMomentumX,
                                        fParticleMomentumY,
                                        fParticleMomentumZ,
                                        0);


        // Generate the particle
        const auto vertex = new G4PrimaryVertex(newParticle.x, newParticle.y, newParticle.z, newParticle.t);
        vertex->SetPrimary(new G4PrimaryParticle(newParticle.pdgid, newParticle.px, newParticle.py, newParticle.pz));
        anEvent->AddPrimaryVertex(vertex);

        // Save all generated particles of the current event
        genParticles.push_back(newParticle);
    }

    // Core function 2: GeneratePrimaryVertex()
    // This is used to set generator parameters
    void MuPARMA::SetNewValue(G4UIcommand *command,
                              G4String value)
    {
        if (command == _ui_pathname)
        {
            // initialization
            startPARMA(value);
        }
        else if (command == _ui_offset)
        {
            auto vec = _ui_offset->GetNew3VectorValue(value);
            fPARMA_additional_config["offset_x"] = vec[0];
            fPARMA_additional_config["offset_y"] = vec[1];
            fPARMA_additional_config["offset_z"] = vec[2];
            resetDimensions();
        }
        else if (command == _ui_shape)
        {
            fPARMA_additional_config["use_shape"] = _ui_shape->GetNewIntValue(value);
            resetDimensions();
        }
        else if (command == _ui_box)
        {
            auto vec = _ui_offset->GetNew3VectorValue(value);
            fPARMA_additional_config["box_lenx"] = vec[0];
            fPARMA_additional_config["box_leny"] = vec[1];
            fPARMA_additional_config["box_lenz"] = vec[2];
            fPARMA_additional_config["offset_z"] = vec[2];
            resetDimensions();
        }
        else if (command == _ui_offset_t_low)
            fPARMA_additional_config["offset_t_low"] = _ui_offset_t_low->GetNewDoubleValue(value);
        else if (command == _ui_offset_t_high)
            fPARMA_additional_config["offset_t_high"] = _ui_offset_t_high->GetNewDoubleValue(value);
        else if (command == _ui_ekin_low)
        {
            fPARMA_additional_config["ekin_cut_low"] = _ui_ekin_low->GetNewDoubleValue(value);
            fPARMAgenerator.emin = static_cast<int>(fPARMA_additional_config["ekin_cut_low"]);
        }
        else if (command == _ui_ekin_high)
        {
            fPARMA_additional_config["ekin_cut_high"] = _ui_ekin_high->GetNewDoubleValue(value);
            fPARMAgenerator.emax = static_cast<int>(fPARMA_additional_config["ekin_cut_high"]);
        }
        else if (command == _ui_particle)
        {
            fPARMA_additional_config["particle_pdgid"] = _ui_particle->GetNewDoubleValue(value);
            fPARMAgenerator.ip = PARMA::pdgid_to_id[static_cast<int>(fPARMA_additional_config["particle_pdgid"])];
        }
        else if (command == _ui_update)
        {
            fPARMAgenerator.UpdateParameters();
        }        
    }

    std::map<std::string, std::string> MuPARMA::getMetaData()
    {
        return {};
    }

    float MuPARMA::getEventWeight()
    {
        return event_weight;
    }
}