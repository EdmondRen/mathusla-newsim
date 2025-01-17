#include <filesystem>
#include <cstdio>
#include <iostream>

#include "MuAnalysis.hh"
#include "G4RunManager.hh"

#include "util.hh"
#include "MuEventAction.hh"
#include "MuPrimaryGeneratorAction.hh"
#include "MuDetectorConstruction.hh"
#include "MuSteppingAction.hh"

namespace Analysis
{
    // Geant4 generated tuple ID
    int tuple_id = 0;

    // Geant4 hit allocator
    G4Allocator<uHit> *HitAllocator = new G4Allocator<uHit>;

    //  ---------------------------------------------------------------------------
    // Default hit class
    uHit::uHit() : G4VHit() {}

    uHit::uHit(G4Step *step, G4TouchableHistory *touchable) : G4VHit()
    {
        (void)touchable;
        const auto track = step->GetTrack();
        const auto step_point = step->GetPostStepPoint();

        this->_particle = track->GetParticleDefinition();
        this->_trackID = track->GetTrackID();
        this->_trackPDG = this->_particle->GetPDGEncoding();
        this->_parentID = track->GetParentID();
        this->_parentPDG = 0; // Fix this later
        this->_edeposit = step->GetTotalEnergyDeposit();
        this->_position = G4LorentzVector(step_point->GetGlobalTime(), step_point->GetPosition());
        this->_momentum = G4LorentzVector(step_point->GetTotalEnergy(), step_point->GetMomentum());

        // G4cout<<"-----hit_x----- "<< _position.x() <<G4endl;

        // Get touchable
        const auto preStepPoint = step->GetPreStepPoint();
        const auto touchable_pre = preStepPoint->GetTouchable();
        // Get copy number for each depth
        G4int totalDepth = touchable_pre->GetHistoryDepth();
        for (int i = 0; i < totalDepth; i++)
            this->_copyNumber.push_back(touchable_pre->GetCopyNumber(i));

        G4ThreeVector worldPos = preStepPoint->GetPosition();
        // Get the transformation matrix from the world to the local coordinate system of the bottom-level touchable
        // and perform the coordinate transform 
        this->_local_coord = touchable_pre->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    }

    //  ---------------------------------------------------------------------------
    // Default detector class
    DefaultDetector::DefaultDetector() : G4VSensitiveDetector("mathusla")
    {
        // Setup the data container
        // clang-format off
        fdata = new util::py::Dict();
        fdata->Add("UID",               util::py::__single__, util::py::__double__); // Unique ID for each event, defined as (run_number*1e9 + event_number). This means each run can have at most 1e9 events.
        fdata->Add("Seed_0",            util::py::__single__, util::py::__int__);
        fdata->Add("Seed_1",            util::py::__single__, util::py::__int__);
        fdata->Add("Hit_x",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_y",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_z",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_t",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_edep",          util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_px",            util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_py",            util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_pz",            util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_trackID",       util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_trackIDparent", util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_pdgID",         util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_pdgIDparent",   util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_isprimary",     util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_processID",     util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_detectorID",    util::py::__vector__, util::py::__double__); // int32 is not long enough. Use double. 
        fdata->Add("Gen_x",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_y",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_z",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_t",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_px",            util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_py",            util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_pz",            util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_pdgID",         util::py::__vector__, util::py::__double__);
        fdata->Add("Gen_index",         util::py::__vector__, util::py::__double__);
        // Step data, optional
        fdata->Add("Step_x",             util::py::__vector__, util::py::__float__);
        fdata->Add("Step_y",             util::py::__vector__, util::py::__float__);
        fdata->Add("Step_z",             util::py::__vector__, util::py::__float__);
        fdata->Add("Step_t",             util::py::__vector__, util::py::__float__);
        fdata->Add("Step_edep",          util::py::__vector__, util::py::__float__);
        fdata->Add("Step_px",            util::py::__vector__, util::py::__float__);
        fdata->Add("Step_py",            util::py::__vector__, util::py::__float__);
        fdata->Add("Step_pz",            util::py::__vector__, util::py::__float__);
        fdata->Add("Step_trackID",       util::py::__vector__, util::py::__int__);
        fdata->Add("Step_trackIDparent", util::py::__vector__, util::py::__int__);
        fdata->Add("Step_pdgID",         util::py::__vector__, util::py::__int__);
        fdata->Add("Step_status",        util::py::__vector__, util::py::__int__);
        // clang-format on
    }

    util::py::Dict *DefaultDetector::GetDataDict()
    {
        return fdata;
    }

    void DefaultDetector::Initialize(G4HCofThisEvent *)
    {
        // Define a hits collection name
        G4String hcName = SensitiveDetectorName + "HitsCollection";
        // Create a hits collection object
        fHitsCollection = new HitsCollection(SensitiveDetectorName, hcName);

        // Get a pointer to runactions
        fMuRunAction = static_cast<const MuRunAction *>(G4RunManager::GetRunManager()->GetUserRunAction());
        this->run_number = fMuRunAction->GetRunNumber();
    }

    G4bool DefaultDetector::ProcessHits(G4Step *step, G4TouchableHistory *touchable)
    {
        // Create a hit
        auto newHit = new uHit(step, touchable);
        // Add the hit in the SD hits collection
        fHitsCollection->insert(newHit);

        return true;
    }

    void DefaultDetector::EndOfEvent(G4HCofThisEvent *)
    {
        const G4Event *event = G4RunManager::GetRunManager()->GetCurrentEvent();
        if (!event)
            return;

        // Get pointer to data dict
        auto &data = *fdata;
        data.clear();

        // Get the per-Event user information
        auto *eventInfo = dynamic_cast<MyEventInformation *>(event->GetUserInformation());
        std::vector<unsigned long> seedInfo = eventInfo->GetInfo();

        // Get the event index
        const auto event_id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        // Get detector construction
        auto detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        auto detectorConstruction_this = dynamic_cast<const MuDetectorConstruction *>(detectorConstruction);
        auto detectorBuilder = detectorConstruction_this->GetBuilder();

        // Don't write down anything if there are less than 4 hits.
        if (fHitsCollection->GetSize() < 4)
            return;

        // Set single values
        data["UID"] = run_number * 1e9 + event_id;
        data["Seed_0"] = *reinterpret_cast<int *>(&seedInfo[1]); // Cast the address of the unsigned long to an int pointer
        data["Seed_1"] = *reinterpret_cast<int *>(&seedInfo[2]); // Cast the address of the unsigned long to an int pointer

        // Process hit collection
        for (std::size_t i = 0; i < fHitsCollection->GetSize(); ++i)
        {
            const auto hit = dynamic_cast<uHit *>(fHitsCollection->GetHit(i));
            data["Hit_x"].push_back(hit->_position.x());
            data["Hit_y"].push_back(hit->_position.y());
            data["Hit_z"].push_back(hit->_position.z());
            data["Hit_t"].push_back(hit->_position.t());
            data["Hit_edep"].push_back(hit->_edeposit);
            data["Hit_px"].push_back(hit->_momentum.px());
            data["Hit_py"].push_back(hit->_momentum.py());
            data["Hit_pz"].push_back(hit->_momentum.pz());
            data["Hit_trackID"].push_back(hit->_trackID);
            data["Hit_trackIDparent"].push_back(hit->_parentID);
            data["Hit_pdgID"].push_back(hit->_trackPDG);
            data["Hit_pdgIDparent"].push_back(hit->_parentPDG);
            data["Hit_processID"].push_back(0); // Fix this later
            data["Hit_detectorID"].push_back(detectorBuilder->GetDetectorID(hit->_copyNumber, hit->_local_coord));

            // Copy number is a list for each hit. Flatten it and separate by -1
            // Start with lowest depth (for example, bar number) to highes depth (detector number)
            // for (int cn : hit->_copyNumber)
            //     data["Hit_copyNumber"].push_back(cn);
            // data["Hit_copyNumber"].push_back(-1); // -1 is used to separate multiple hits in one event.
        }

        // Process generated particles
        auto GenPrticles = MuPrimaryGeneratorAction::GetLastEvent();
        for (auto particle : GenPrticles)
        {
            data["Gen_x"].push_back(particle.x);
            data["Gen_y"].push_back(particle.y);
            data["Gen_z"].push_back(particle.z);
            data["Gen_t"].push_back(particle.t);
            data["Gen_px"].push_back(particle.px);
            data["Gen_py"].push_back(particle.py);
            data["Gen_pz"].push_back(particle.pz);
            data["Gen_pdgID"].push_back(particle.pdgid);
            data["Gen_index"].push_back(particle.index);
        }

        // Process Step data (only when saving step is enabled)
        auto userStepAction = dynamic_cast<const MuSteppingAction *>(G4RunManager::GetRunManager()->GetUserSteppingAction());
        if (userStepAction->ENABLE_STEPS)
        {
            auto fStepDataStore = userStepAction->fStepDataStore;
            for (size_t ii = 0; ii < fStepDataStore->_step_x.size(); ii++)
            {
                data["Step_x"].push_back(fStepDataStore->_step_x[ii]);
                data["Step_y"].push_back(fStepDataStore->_step_y[ii]);
                data["Step_z"].push_back(fStepDataStore->_step_z[ii]);
                data["Step_t"].push_back(fStepDataStore->_step_t[ii]);
                data["Step_edep"].push_back(fStepDataStore->_step_edep[ii]);
                data["Step_px"].push_back(fStepDataStore->_step_px[ii]);
                data["Step_py"].push_back(fStepDataStore->_step_py[ii]);
                data["Step_pz"].push_back(fStepDataStore->_step_pz[ii]);
                data["Step_trackID"].push_back(fStepDataStore->_step_trackID[ii]);
                data["Step_trackIDparent"].push_back(fStepDataStore->_step_trackIDparent[ii]);
                data["Step_pdgID"].push_back(fStepDataStore->_step_pdg[ii]);
                data["Step_status"].push_back(fStepDataStore->_step_status[ii]);
            }
        }

        // Fill them into the tuple
        FillNTuple(data);
    }

    //  ---------------------------------------------------------------------------
    // Common analysis helper functions
    // Setup ROOT Analysis Tool
    void Setup()
    {
        // #if (G4VERSION_NUMBER <= 1100)
        //         delete G4AnalysisManager::Instance();
        // #else
        //         G4AnalysisManager::Instance()->Clear();
        // #endif
    }

    // Open Output File
    bool Open(const std::string &path)
    {
        auto path_full = path + ".root";
        // util::py::print("File to create", path_full);
        // if (std::filesystem::exists(path_full)){
        //     std::remove(path_full.c_str());
        //     util::py::print("File alreay exists, removed.", path_full);
        // }
        return G4AnalysisManager::Instance()->OpenFile(path_full);
    }

    // Save Output
    bool Save()
    {
        return G4AnalysisManager::Instance()->Write() && G4AnalysisManager::Instance()->CloseFile();
    }

    // Create ROOT NTuple
    bool MuCreateNTuple(util::py::Dict &data, const std::string &name)
    {
        const auto manager = G4AnalysisManager::Instance();
        tuple_id = manager->CreateNtuple(name, name);
        const auto size = data.size();

        for (int i = 0; i < size; i++)
        {
            auto &data_column = data[i];
            switch (data_column.key_type)
            {
            case util::py::__single__:
                switch (data_column.data_type)
                {
                case util::py::__int__:
                    manager->CreateNtupleIColumn(tuple_id, data_column.key_name);
                    break;
                case util::py::__float__:
                    manager->CreateNtupleFColumn(tuple_id, data_column.key_name);
                    break;
                case util::py::__double__:
                    manager->CreateNtupleDColumn(tuple_id, data_column.key_name);
                    break;
                case util::py::__string__:
                    manager->CreateNtupleSColumn(tuple_id, data_column.key_name);
                    break;
                }
                break;
            case util::py::__vector__:
                switch (data_column.data_type)
                {
                case util::py::__int__:
                    manager->CreateNtupleIColumn(tuple_id, data_column.key_name, data_column.data_vec_int);
                    break;
                case util::py::__float__:
                    manager->CreateNtupleFColumn(tuple_id, data_column.key_name, data_column.data_vec_float);
                    break;
                case util::py::__double__:
                    manager->CreateNtupleDColumn(tuple_id, data_column.key_name, data_column.data_vec_double);
                    break;
                default:
                    break;
                }
                break;
            }
        }

        manager->FinishNtuple(tuple_id);
        return true;
    }

    // Fill ROOT NTuple
    bool FillNTuple(util::py::Dict &data)
    {
        const auto manager = G4AnalysisManager::Instance();
        const auto size = data.size();

        for (int i = 0; i < size; i++)
        {
            auto &data_column = data[i];
            switch (data_column.key_type)
            {
            case util::py::__single__:
                switch (data_column.data_type)
                {
                case util::py::__int__:
                    manager->FillNtupleIColumn(tuple_id, data_column.key_index, data_column.data_int);
                    break;
                case util::py::__float__:
                    manager->FillNtupleFColumn(tuple_id, data_column.key_index, data_column.data_float);
                    break;
                case util::py::__double__:
                    manager->FillNtupleDColumn(tuple_id, data_column.key_index, data_column.data_double);
                    break;
                case util::py::__string__:
                    manager->FillNtupleSColumn(tuple_id, data_column.key_index, data_column.data_string);
                    break;
                }
                break;
            default:
                break;
            }
        }

        manager->AddNtupleRow(tuple_id);
        return true;
    }

} // namespace Analysis