
#include "MuAnalysis.hh"

namespace Analysis
{

    //__Setup ROOT Analysis Tool____________________________________________________________________
    void Setup()
    {
        // #if (G4VERSION_NUMBER <= 1100)
        //         delete G4AnalysisManager::Instance();
        // #else
        //         G4AnalysisManager::Instance()->Clear();
        // #endif
        auto MuAnalysisManager = G4AnalysisManager::Instance();

        MuAnalysisManager->SetNtupleMerging(false);
        MuAnalysisManager->SetVerboseLevel(0);

        

    }

    //__Open Output File____________________________________________________________________________
    bool Open(const std::string &path)
    {
        return G4AnalysisManager::Instance()->OpenFile(path);
    }

    //__Save Output_________________________________________________________________________________
    bool Save()
    {
        return G4AnalysisManager::Instance()->Write() && G4AnalysisManager::Instance()->CloseFile();
    }

    //__Create ROOT NTuple__________________________________________________________________________
    bool CreateNTuple(KEY_DATA_t &data, const std::string &name)
    {
        const auto manager = G4AnalysisManager::Instance();
        const auto id = manager->CreateNtuple(name, name);
        const auto size = data.size();

        for (const auto& data_column : data)
        {   
            switch(data_column.key_type){
                case SINGLE:
                    switch (data_column.data_type)
                    {
                    case TUPLE_INT:
                        manager->CreateNtupleDColumn(id, data_column.key_name);
                        break;
                    case TUPLE_FLOAT:
                        manager->CreateNtupleIColumn(id, data_column.key_name);
                        break;
                    case TUPLE_DOUBLE:
                        manager->CreateNtupleFColumn(id, data_column.key_name);
                        break;
                    case TUPLE_STRING:
                        manager->CreateNtupleSColumn(id, data_column.key_name);
                        break;                                                                        
                    }
                    break;
                case VECTOR:
                    switch (data_column.data_type)
                    {
                    case TUPLE_INT:
                        manager->CreateNtupleDColumn(id, data_column.key_name, *data_column.vec_int);
                        break;
                    case TUPLE_FLOAT:
                        manager->CreateNtupleIColumn(id, data_column.key_name);
                        break;
                    case TUPLE_DOUBLE:
                        manager->CreateNtupleFColumn(id, data_column.key_name);
                        break;
                    }
                    break;                
            }
        }
                    

        manager->FinishNtuple(id);
        return _ntuple.insert({name, id}).second;
    }
    //----------------------------------------------------------------------------------------------

    //__Fill ROOT NTuple____________________________________________________________________________
    bool FillNTuple(const std::string &name,
                    const DataKeyTypeList &types,
                    const DataEntry &single_values,
                    const DataEntryList &vector_values)
    {

        const auto search = _ntuple.find(name);
        if (search == _ntuple.cend())
            return false;

        auto &data = _ntuple_data[name];
        const auto vector_size = vector_values.size();
        if (data.size() != vector_size)
        {
            std::cout << "data size" << data.size() << std::endl;
            std::cout << "vector size" << vector_size << std::endl;
            return false;
        }

        for (std::size_t i{}; i < vector_size; ++i)
            data[i] = vector_values[i];

        const auto id = search->second;
        const auto manager = G4AnalysisManager::Instance();
        const auto size = types.size();
        for (std::size_t index{}, single_index{}; index < size; ++index)
        {
            if (types[index] == DataKeyType::Single)
            {
                manager->FillNtupleDColumn(id, index, single_values[single_index++]);
            }
        }

        manager->AddNtupleRow(id);
        return true;
    }
} // namespace Analysis