// run with: root -l -q CreateSampleTree.C

void CreateSampleTree() {
    TFile file("source.root", "RECREATE");
    TTree tree("sourceTree", "Sample Tree");

    Int_t x;
    Double_t y;
    tree.Branch("x", &x);
    tree.Branch("y", &y);

    for (Int_t i = 0; i < 100; i++) {
        x = i;
        y = i * 0.5;
        tree.Fill();
    }

    tree.Write();
    file.Close();
}