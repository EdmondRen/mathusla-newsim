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

void mergeTrees(std::string source_fname, std::string target_fname, std::string tree_name)
{
    // Open the source and target files
    TFile *source_file = TFile::Open(source_fname.c_str(), "READ");
    TFile *target_file = TFile::Open(target_fname.c_str(), "UPDATE");

    // If tree name not given, auto get Tree name
    source_file->GetListOfKeys();
    TKey *keyref;
    TIter nextref(source_file->GetListOfKeys());
    if (tree_name.size() == 0)
    {
        while ((keyref = (TKey *)nextref()))
        {
            if (TString(keyref->GetClassName()).BeginsWith("TTree"))
            {
                tree_name = keyref->GetName();
                std::cout << "Tree name not provided. Merging the first tree with name: " << tree_name << std::endl;
                break;
            }
        }
    }

    // Make a counter for number of files merged
    TParameter<int> *N_MERGED = nullptr;
    if (source_file->GetListOfKeys()->FindObject("N_MERGED"))
        N_MERGED = (TParameter<int> *)source_file->Get("N_MERGED");
    if (!N_MERGED)
        N_MERGED = new TParameter<int>("N_MERGED", 1);

    N_MERGED->SetVal(N_MERGED->GetVal() + 1);

    // Retrieve the trees
    TTree *source_tree = (TTree *)source_file->Get(tree_name.c_str());
    TTree *target_tree = (TTree *)target_file->Get(tree_name.c_str());

    // Method 1: use TChain
    // // Add trees to a list
    // TList *list = new TList;
    // list->Add(source_tree);
    // list->Add(target_tree);

    // // Merge the trees
    // if (source_tree && target_tree)
    //     target_tree->Merge(list);
    // else
    //     std::cerr << "Error: One or both trees not found!" << std::endl;

    // Method 2: direct copy
    // Clone the tree structure (but not the entries)
    target_tree->CopyAddresses(source_tree);
    target_tree->CopyEntries(source_tree);
    // // Loop over entries in the external tree and copy them to the current tree
    // Long64_t nEntries = source_tree->GetEntries();
    // for (Long64_t i = 0; i < nEntries; i++) {
    //     source_tree->GetEntry(i);  // Load entry data
    //     target_tree->Fill();        // Fill current tree with loaded entry
    // }

    // Write and close the files
    target_file->cd();
    // Write the Counter to the file
    N_MERGED->Write("", TObject::kOverwrite);
    target_file->Write("", TObject::kOverwrite);
    target_file->Close();
    source_file->Close();
}

int main(int argc, char *argv[])
{   

    // Handel input arguments
    std::string source_fname, target_fname, tree_name;
    if (argc < 3 || argc > 5)
        throw std::invalid_argument("Usage: ./merge source target [tree_name]");

    else if (argc >= 3)
    {
        source_fname = argv[1];
        target_fname = argv[2];
    }

    if (argc == 4)
        tree_name = argv[3];

    std::cout << "Source: " << source_fname << std::endl;
    std::cout << "Target: " << target_fname << std::endl;
    std::cout << "Tree name: " << tree_name << std::endl;


    if (std::filesystem::exists(target_fname)) {
        mergeTrees(source_fname, target_fname, tree_name);
        std::cout << "File merging completed" << std::endl;
    } else {
        std::filesystem::copy_file(source_fname, target_fname, std::filesystem::copy_options::overwrite_existing);
        std::cout << "Target file does not exist. Source file copied as target.\n";
    }


    return 0;
}