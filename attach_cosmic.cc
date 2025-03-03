#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TKey.h>
#include <TTree.h>
#include <TChain.h>
#include <TParameter.h>
#include <TRandom3.h>

TRandom3 rng;

void AttachTrees(std::string signal_fname, std::string cosmic_fname, std::string output_name, float poisson_mean)
{

    rng.SetSeed();

    // Open the source and target files
    std::string tree_name = "data";
    TFile *signal_file = TFile::Open(signal_fname.c_str(), "READ");
    TFile *cosmic_file = TFile::Open(cosmic_fname.c_str(), "READ");
    TFile *out_file = TFile::Open(output_name.c_str(), "RECREATE");

    // Retrieve the trees
    TTree *t1 = (TTree *)signal_file->Get(tree_name.c_str());
    TTree *t2 = (TTree *)cosmic_file->Get(tree_name.c_str());
    TTree *outTree = t1->CloneTree(0); // Clone structure but not the entries
    // Don't forget the metadata
    TTree *treeMetadata = (TTree *)signal_file->Get("metadata");
    treeMetadata->CloneTree()->Write("metadata", TObject::kOverwrite);

    // Create maps for storing branch addresses
    std::vector<TBranch *> branches1, branches2;
    // std::vector<void*> addresses1, addresses2;
    static const int N = 1024;
    void *addresses1[N];
    void *addresses2[N];
    std::vector<float> *branchData_vfloat[N];
    std::vector<double> *branchData_vdouble[N];
    std::vector<int> *branchData_vint[N];
    std::vector<Long64_t> *branchData_vint64[N];    

    TObjArray *branchList = t1->GetListOfBranches();
    for (int i = 0; i < branchList->GetEntries(); ++i)
    {
        TBranch *br1 = (TBranch *)branchList->At(i);
        TBranch *br2 = t2->GetBranch(br1->GetName());
        if (!br2)
        {
            std::cerr << "Error: Branch " << br1->GetName() << " not found in second tree." << std::endl;
            return;
        }

        branches1.push_back(br1);
        branches2.push_back(br2);

        // Determine branch type
        void *addr1 = nullptr;
        void *addr2 = nullptr;

        TClass *cls;
        EDataType type;
        br1->GetExpectedType(cls, type);

        if (cls)
        { // Assume it's a vector if it's an object
            if (cls->GetName() == std::string("vector<int>"))
            {
                auto _addr1 = new std::vector<int>();
                auto _addr2 = new std::vector<int>();
                branchData_vint[i] = _addr1;
                branchData_vint[N - i] = _addr2;
                t1->SetBranchAddress(br1->GetName(), &branchData_vint[i]);
                t2->SetBranchAddress(br2->GetName(), &branchData_vint[N - i]);
                addr1 = (void *)_addr1;
                addr2 = (void *)_addr2;
            }
            else if (cls->GetName() == std::string("vector<float>"))
            {
                auto _addr1 = new std::vector<float>();
                auto _addr2 = new std::vector<float>();
                branchData_vfloat[i] = _addr1;
                branchData_vfloat[N - i] = _addr2;
                t1->SetBranchAddress(br1->GetName(), &branchData_vfloat[i]);
                t2->SetBranchAddress(br2->GetName(), &branchData_vfloat[N - i]);
                addr1 = (void *)_addr1;
                addr2 = (void *)_addr2;
            }
            else if (cls->GetName() == std::string("vector<double>"))
            {
                auto _addr1 = new std::vector<double>();
                auto _addr2 = new std::vector<double>();
                branchData_vdouble[i] = _addr1;
                branchData_vdouble[N - i] = _addr2;
                t1->SetBranchAddress(br1->GetName(), &branchData_vdouble[i]);
                t2->SetBranchAddress(br2->GetName(), &branchData_vdouble[N - i]);
                addr1 = (void *)_addr1;
                addr2 = (void *)_addr2;
            }
            else if (cls->GetName() == std::string("vector<Long64_t>"))
            {
                auto _addr1 = new std::vector<Long64_t>();
                auto _addr2 = new std::vector<Long64_t>();
                branchData_vint64[i] = _addr1;
                branchData_vint64[N - i] = _addr2;
                t1->SetBranchAddress(br1->GetName(), &branchData_vint64[i]);
                t2->SetBranchAddress(br2->GetName(), &branchData_vint64[N - i]);
                addr1 = (void *)_addr1;
                addr2 = (void *)_addr2;
            }
            else
            {
                std::cerr << "Error: Unsupported vector type for branch " << br1->GetName() << ", type: " << cls->GetName() << std::endl;
                return;
            }
            std::cerr << "Branch " << br1->GetName() << ", type: " << cls->GetName() << std::endl;
        }
        else
        { // Assume primitive type
            addr1 = operator new(br1->GetTotalSize());
            addr2 = operator new(br2->GetTotalSize());
            std::cerr << "Branch " << br1->GetName() << std::endl;

            t1->SetBranchAddress(br1->GetName(), addr1);
            t2->SetBranchAddress(br2->GetName(), addr2);
        }

        addresses1[i] = addr1;
        addresses2[i] = addr2;

        // Don't do this. CloneTree() already set the address of outTree.
        // outTree->SetBranchAddress(br1->GetName(), addr1);
    }

    // Append the content of the second file to the first one.
    Long64_t nEntries = t1->GetEntries();
    Long64_t entry_cosmic = 0;
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t1->GetEntry(i);

        // Sample to see if want to add the cosmic ray to this event or not
        int k = rng.Poisson(poisson_mean);
        for (;k>0;k--, entry_cosmic++)
        {

            t2->GetEntry(entry_cosmic);

            for (size_t j = 0; j < branches1.size(); ++j)
            {
                TClass *cls;
                EDataType type;
                branches1[j]->GetExpectedType(cls, type);

                if (cls)
                { // If it's a vector, concatenate
                    if (cls->GetName() == std::string("vector<int>"))
                    {
                        auto *v1 = static_cast<std::vector<int> *>(addresses1[j]);
                        auto *v2 = static_cast<std::vector<int> *>(addresses2[j]);
                        v1->insert(v1->end(), v2->begin(), v2->end());
                    }
                    else if (cls->GetName() == std::string("vector<float>"))
                    {
                        auto *v1 = static_cast<std::vector<float> *>(addresses1[j]);
                        auto *v2 = static_cast<std::vector<float> *>(addresses2[j]);
                        v1->insert(v1->end(), v2->begin(), v2->end());
                    }
                    else if (cls->GetName() == std::string("vector<double>"))
                    {
                        auto *v1 = static_cast<std::vector<double> *>(addresses1[j]);
                        auto *v2 = static_cast<std::vector<double> *>(addresses2[j]);
                        v1->insert(v1->end(), v2->begin(), v2->end());
                    }
                    else if (cls->GetName() == std::string("vector<Long64_t>"))
                    {
                        auto *v1 = static_cast<std::vector<Long64_t> *>(addresses1[j]);
                        auto *v2 = static_cast<std::vector<Long64_t> *>(addresses2[j]);
                        v1->insert(v1->end(), v2->begin(), v2->end());
                    }
                }
                else
                { // Copy primitive types directly
                    std::memcpy(addresses1[j], addresses2[j], branches1[j]->GetTotalSize());
                }
            }
        }


        outTree->Fill();
    }
    out_file->Write("", TObject::kOverwrite);
    out_file->Close();
    signal_file->Close();
    cosmic_file->Close();
}

int main(int argc, char *argv[])
{

    // Handel input arguments
    std::string signal_fname, cosmic_fname, output_fname;
    float prob;
    if (argc != 5)
        throw std::invalid_argument("Usage: ./attach signal_file cosmic_file output_file mean_cosmic_per_event");

    signal_fname = argv[1];
    cosmic_fname = argv[2];
    output_fname = argv[3];
    prob = std::stof(argv[4]);

    std::cout << "signal_file: " << signal_fname << std::endl;
    std::cout << "cosmic_file: " << cosmic_fname << std::endl;
    std::cout << "output_file: " << output_fname << std::endl;

    AttachTrees(signal_fname, cosmic_fname, output_fname, prob);

    std::cout << "File merging completed" << std::endl;

    return 0;
}