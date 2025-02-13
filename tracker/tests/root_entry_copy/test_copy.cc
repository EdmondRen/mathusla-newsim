#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <TObjArray.h>

#include <vector>
#include <string>
#include <iostream>
#include <variant>

void CopyTreeBranches(TTree *sourceTree, TTree *destTree)
{
    // Get the list of branches in the source tree
    TObjArray *branches = sourceTree->GetListOfBranches();

    // Create a vector to hold pointers to the variables used for branch data
    std::vector<std::variant<
        std::unique_ptr<Int_t>,
        std::unique_ptr<Float_t>,
        std::unique_ptr<Double_t>,
        std::unique_ptr<Bool_t>,
        std::unique_ptr<Char_t>,
        std::unique_ptr<std::vector<float>>,
        std::unique_ptr<std::vector<double>>,
        std::unique_ptr<std::vector<int>>>>
        branchData;

    std::vector<std::variant<
        std::vector<float> *,
        std::vector<double> *,
        std::vector<int> *>>
        branchData_vec;

    const int n = 200;
    std::vector<float> *branchData_vfloat[n];
    std::vector<double> *branchData_vdouble[n];
    std::vector<int> *branchData_vint[n];
    int vfloat_n = 0;
    int vdouble_n = 0;
    int vint_n = 0;

    // Iterate over all branches in the source tree
    for (Int_t i = 0; i < branches->GetEntries(); i++)
    {
        TBranch *branch = (TBranch *)branches->At(i);
        const char *branchName = branch->GetName();

        // Get the leaf associated with the branch
        TLeaf *leaf = branch->GetLeaf(branchName);
        if (!leaf)
            continue;

        // Determine the type of the leaf
        const char *leafType = leaf->GetTypeName();
        std::cout << "  Leaf type: " << leafType << " for branch " << branchName << std::endl;

        // Create a new branch in the destination tree with the same name and type
        if (strcmp(leafType, "Int_t") == 0)
        {
            auto value = std::make_unique<Int_t>();
            destTree->Branch(branchName, value.get());
            sourceTree->SetBranchAddress(branchName, value.get());
            branchData.push_back(std::move(value));
        }
        else if (strcmp(leafType, "Float_t") == 0)
        {
            auto value = std::make_unique<Float_t>();
            destTree->Branch(branchName, value.get());
            sourceTree->SetBranchAddress(branchName, value.get());
            branchData.push_back(std::move(value));
        }
        else if (strcmp(leafType, "Double_t") == 0)
        {
            auto value = std::make_unique<Double_t>();
            destTree->Branch(branchName, value.get());
            sourceTree->SetBranchAddress(branchName, value.get());
            branchData.push_back(std::move(value));
        }
        else if (strcmp(leafType, "Bool_t") == 0)
        {
            auto value = std::make_unique<Bool_t>();
            destTree->Branch(branchName, value.get());
            sourceTree->SetBranchAddress(branchName, value.get());
            branchData.push_back(std::move(value));
        }
        else if (strcmp(leafType, "Char_t") == 0)
        {
            auto value = std::make_unique<Char_t>();
            destTree->Branch(branchName, value.get());
            sourceTree->SetBranchAddress(branchName, value.get());
            branchData.push_back(std::move(value));
        }

        else if (strcmp(leafType, "vector<float>") == 0)
        {
            auto value = new std::vector<float>();
            branchData_vfloat[vfloat_n] = value;
            destTree->Branch(branchName, &(*branchData_vfloat[vfloat_n]));          // Use the raw pointer
            sourceTree->SetBranchAddress(branchName, &branchData_vfloat[vfloat_n]); // Use the raw pointer
            vfloat_n += 1;

            // Convert raw pointer to unique_ptr and save to the list
            std::unique_ptr<std::vector<float>> uniquePtr(value);
            branchData.push_back(std::move(uniquePtr));
        }
        else if (strcmp(leafType, "vector<double>") == 0)
        {
            auto value = new std::vector<double>();
            branchData_vdouble[vdouble_n] = value;
            destTree->Branch(branchName, &(*branchData_vdouble[vdouble_n]));          // Use the raw pointer
            sourceTree->SetBranchAddress(branchName, &branchData_vdouble[vdouble_n]); // Use the raw pointer
            vdouble_n += 1;

            // Convert raw pointer to unique_ptr and save to the list
            std::unique_ptr<std::vector<double>> uniquePtr(value);
            branchData.push_back(std::move(uniquePtr));            
        }
        else if (strcmp(leafType, "vector<int>") == 0)
        {
            auto value = new std::vector<int>();
            branchData_vint[vint_n] = value;
            destTree->Branch(branchName, &(*branchData_vint[vint_n]));          // Use the raw pointer
            sourceTree->SetBranchAddress(branchName, &branchData_vint[vint_n]); // Use the raw pointer
            vint_n += 1;

            // Convert raw pointer to unique_ptr and save to the list
            std::unique_ptr<std::vector<int>> uniquePtr(value);
            branchData.push_back(std::move(uniquePtr));            
        }

        else
        {
            // Handle other types or arrays if necessary
            std::cerr << "Unsupported type: " << leafType << " for branch " << branchName << std::endl;
            continue;
        }
    }

    // Copy the data from the source tree to the destination tree
    Long64_t nEntries = sourceTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++)
    {
        // std::cout << "entry " << i << std::endl;
        sourceTree->GetEntry(i);
        destTree->Fill();
    }
}

int main()
{
    // Open the source file and get the source tree
    TFile *sourceFile = TFile::Open("../../../build/data/run_0.root");
    TTree *sourceTree = (TTree *)sourceFile->Get("data");

    // Create a new file and a new tree
    TFile *destFile = new TFile("../../../build/data/destination.root", "RECREATE");
    TTree *destTree = new TTree("data", "Tree with copied branches");

    // Copy branches from sourceTree to destTree
    CopyTreeBranches(sourceTree, destTree);
    // CopyBranches(sourceTree, destTree);

    // Write the destination tree to the file
    destFile->Write();
    destFile->Close();

    // Close the source file
    sourceFile->Close();

    return 0;
}