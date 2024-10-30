#include <filesystem>
#include <cstdio>
#include <iostream>


#include "MuAnalysis.hh"
#include "G4RunManager.hh"

#include "util.hh"

namespace Analysis
{
    int tuple_id = 0;

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
        auto path_full = path+".root";
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
    bool CreateNTuple(util::py::Dict &data, const std::string &name)
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

    // Geant4 hit allocator
    G4Allocator<uHit> *HitAllocator = new G4Allocator<uHit>;

    //  ---------------------------------------------------------------------------
    // Default hit class
    uHit::uHit() : G4VHit()
    {
    }

    uHit::uHit(G4Step *step) : G4VHit()
    {
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
    }

    //  ---------------------------------------------------------------------------
    // Default detector class
    DefaultDetector::DefaultDetector() : G4VSensitiveDetector("mathusla")
    {
        // Setup the data container
        // clang-format off
        fdata = new util::py::Dict();
        fdata->Add("Entry_generated",   util::py::__single__, util::py::__float__);
        fdata->Add("Hit_x",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_y",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_z",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_t",             util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_edep",          util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_px",            util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_py",            util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_pz",            util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_trackID",       util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_trackIDparent", util::py::__vector__, util::py::__float__);
        fdata->Add("Hit_pdgID",         util::py::__vector__, util::py::__double__);
        fdata->Add("Hit_pdgIDparent",   util::py::__vector__, util::py::__double__);
        fdata->Add("Hit_isprimary",     util::py::__vector__, util::py::__int__);
        fdata->Add("Hit_processID",     util::py::__vector__, util::py::__int__);
        fdata->Add("Gen_x",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_y",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_z",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_t",             util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_px",            util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_py",            util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_pz",            util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_edep",          util::py::__vector__, util::py::__float__);
        fdata->Add("Gen_pdgID",         util::py::__vector__, util::py::__double__);
        // clang-format on
    }

    void DefaultDetector::Initialize(G4HCofThisEvent *event)
    {
        // Define a hits collection name
        G4String hcName = SensitiveDetectorName + "HitsCollection";
        // Create a hits collection object
        fHitsCollection = new HitsCollection(SensitiveDetectorName, hcName);
    }

    G4bool DefaultDetector::ProcessHits(G4Step *step, G4TouchableHistory *touchable)
    {
        // Create a hit
        auto newHit = new uHit(step);
        // Add the hit in the SD hits collection
        fHitsCollection->insert(newHit);
    }

    void DefaultDetector::EndOfEvent(G4HCofThisEvent *)
    {
        // // Get the event index
        // const auto event_id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        // // Process hit collection
        // auto &data = *fdata;
        // data.clear();
        // for (std::size_t i = 0; i < fHitsCollection->GetSize(); ++i)
        // {
        //     const auto hit = dynamic_cast<uHit *>(fHitsCollection->GetHit(i));
        //     data["Hit_x"].push_back(hit->_position.x());
        //     data["Hit_y"].push_back(hit->_position.y());
        //     data["Hit_z"].push_back(hit->_position.z());
        //     data["Hit_t"].push_back(hit->_position.t());
        //     data["Hit_edep"].push_back(hit->_edeposit);
        //     data["Hit_px"].push_back(hit->_momentum.px());
        //     data["Hit_py"].push_back(hit->_momentum.py());
        //     data["Hit_pz"].push_back(hit->_momentum.pz());
        //     data["Hit_trackID"].push_back(hit->_trackID);            
        //     data["Hit_trackIDparent"].push_back(hit->_parentID);            
        //     data["Hit_pdgID"].push_back(hit->_trackPDG);            
        //     data["Hit_pdgIDparent"].push_back(hit->_parentPDG);            
        //     data["Hit_processID"].push_back(0); //Fix this later                
        // }

        // // Fill them into the tuple
        // FillNTuple(data);
    }

} // namespace Analysis