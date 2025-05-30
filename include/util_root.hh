// myheader.h
#ifndef util_root_hh
#define util_root_hh

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <variant>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TStyle.h>

// ROOT plot related
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TColor.h>

namespace iroot
{
    namespace style
    {
        static void setTableauPalette()
        {
            // Define RGB values for Tableau colors
            const int nColors = 10;
            double tableauColors[nColors][3] = {
                {31 / 255.0, 119 / 255.0, 180 / 255.0},  // Blue
                {255 / 255.0, 127 / 255.0, 14 / 255.0},  // Orange
                {44 / 255.0, 160 / 255.0, 44 / 255.0},   // Green
                {214 / 255.0, 39 / 255.0, 40 / 255.0},   // Red
                {148 / 255.0, 103 / 255.0, 189 / 255.0}, // Purple
                {140 / 255.0, 86 / 255.0, 75 / 255.0},   // Brown
                {227 / 255.0, 119 / 255.0, 194 / 255.0}, // Pink
                {127 / 255.0, 127 / 255.0, 127 / 255.0}, // Gray
                {188 / 255.0, 189 / 255.0, 34 / 255.0},  // Olive
                {23 / 255.0, 190 / 255.0, 207 / 255.0}   // Cyan
            };

            // Create ROOT colors
            std::vector<int> colorIndices;
            for (int i = 0; i < nColors; ++i)
            {
                int colorIndex = TColor::GetFreeColorIndex();
                new TColor(colorIndex, tableauColors[i][0], tableauColors[i][1], tableauColors[i][2]);
                colorIndices.push_back(colorIndex);
            }

            // Set the palette
            gStyle->SetPalette(colorIndices.size(), &colorIndices[0]);
        }

        inline void setStyle()
        {
            gROOT->SetStyle("Modern");
            setTableauPalette();
        }
    } // namespace style

    /*
    namespace plt
    {

        class figure
        {
        private:
            std::unique_ptr<TCanvas> canvas;                                                                                     // The canvas for drawing
            int colorIndex;                                                                                                      // Current color index for automatic color management
            std::vector<std::variant<std::unique_ptr<TGraph>, std::unique_ptr<TGraphAsymmErrors>, std::unique_ptr<TH1D>>> plots; // Container for plot objects

        public:
            // Constructor: creates a canvas and initializes the color index
            figure(const std::string &title = "Figure", int width = 800, int height = 600)
                : canvas(std::make_unique<TCanvas>("canvas", title.c_str(), width, height)), colorIndex(0)
            {
                gStyle->SetPalette(kRainBow); // Use a nice color palette
            }

            // Set axis limits
            void xlim(double xmin, double xmax)
            {
                canvas->cd();
                canvas->GetXaxis()->SetLimits(xmin, xmax);
            }

            void ylim(double ymin, double ymax)
            {
                canvas->cd();
                canvas->GetYaxis()->SetRangeUser(ymin, ymax);
            }

            // Set axis labels
            void xlabel(const std::string &label)
            {
                canvas->cd();
                canvas->GetXaxis()->SetTitle(label.c_str());
            }

            void ylabel(const std::string &label)
            {
                canvas->cd();
                canvas->GetYaxis()->SetTitle(label.c_str());
            }

            // Set figure title
            void title(const std::string &title)
            {
                canvas->cd();
                canvas->SetTitle(title.c_str());
            }

            // Set axis scale (linear or log)
            void xscale(const std::string &scale)
            {
                canvas->cd();
                if (scale == "log")
                    canvas->SetLogx(1);
                else if (scale == "linear")
                    canvas->SetLogx(0);
                else
                    throw std::invalid_argument("Invalid scale type. Use 'log' or 'linear'.");
            }

            void yscale(const std::string &scale)
            {
                canvas->cd();
                if (scale == "log")
                    canvas->SetLogy(1);
                else if (scale == "linear")
                    canvas->SetLogy(0);
                else
                    throw std::invalid_argument("Invalid scale type. Use 'log' or 'linear'.");
            }

            // Plot x-y data using TGraph
            void plot(const std::vector<double> &x, const std::vector<double> &y, const std::string &style = "l", const std::string &label = "", int color = -1)
            {
                if (x.size() != y.size())
                    throw std::invalid_argument("x and y vectors must have the same size.");

                canvas->cd();
                auto graph = std::make_unique<TGraph>(x.size(), x.data(), y.data());

                // Set style and color
                if (color == -1)
                    color = GetNextColor();
                graph->SetLineColor(color);
                graph->SetMarkerColor(color);

                if (style.find("l") != std::string::npos)
                    graph->SetLineWidth(2);
                if (style.find("p") != std::string::npos)
                    graph->SetMarkerStyle(20);

                // Add label if provided
                if (!label.empty())
                    graph->SetTitle(label.c_str());

                // Draw the graph
                graph->Draw((style + " SAME").c_str());

                // Store the graph in the container
                plots.push_back(std::move(graph));
                canvas->Update();
            }

            // Plot x-y data with errors using TGraphAsymmErrors
            void errorbar(const std::vector<double> &x, const std::vector<double> &y,
                          const std::vector<double> &xerr, const std::vector<double> &yerr,
                          const std::string &style = "p", const std::string &label = "", int color = -1)
            {
                if (x.size() != y.size() || x.size() != xerr.size() || x.size() != yerr.size())
                    throw std::invalid_argument("x, y, xerr, and yerr vectors must have the same size.");

                canvas->cd();
                auto graph = std::make_unique<TGraphAsymmErrors>(x.size(), x.data(), y.data(), xerr.data(), xerr.data(), yerr.data(), yerr.data());

                // Set style and color
                if (color == -1)
                    color = GetNextColor();
                graph->SetLineColor(color);
                graph->SetMarkerColor(color);

                if (style.find("p") != std::string::npos)
                    graph->SetMarkerStyle(20);

                // Add label if provided
                if (!label.empty())
                    graph->SetTitle(label.c_str());

                // Draw the graph
                graph->Draw("APZ SAME");

                // Store the graph in the container
                plots.push_back(std::move(graph));
                canvas->Update();
            }

            // Create and plot a histogram using TH1D
            void hist(const std::vector<double> &x, int bins = 10, const std::string &style = "hist", const std::string &label = "", int color = -1, double xmin = 0, double xmax = 0)
            {
                if (x.empty())
                    throw std::invalid_argument("Input data vector x is empty.");

                canvas->cd();

                // Determine range
                if (xmin == xmax)
                    xmin = *std::min_element(x.begin(), x.end());
                    xmax = *std::max_element(x.begin(), x.end());

                auto h = std::make_unique<TH1D>("h", label.c_str(), bins, xmin, xmax);

                // Fill the histogram
                for (double val : x)
                    h->Fill(val);

                // Set color
                if (color == -1)
                    color = GetNextColor();
                h->SetLineColor(color);
                h->SetFillColor(color);

                // Set style
                std::string drawOption = style + " SAME";

                // Draw the histogram
                h->Draw(drawOption.c_str());

                // Store the histogram in the container
                plots.push_back(std::move(h));
                canvas->Update();
            }

            // Plot a TH1D as a stair-step plot
            void stair(TH1D *h, const std::string &style = "hist", const std::string &label = "", int color = -1)
            {
                if (!h)
                    throw std::invalid_argument("Histogram pointer is null.");

                canvas->cd();

                // Set color
                if (color == -1)
                    color = GetNextColor();
                h->SetLineColor(color);

                // Add label if provided
                if (!label.empty())
                    h->SetTitle(label.c_str());

                // Set style
                std::string drawOption = style + " SAME";

                // Draw the histogram
                h->Draw(drawOption.c_str());

                canvas->Update();
            }

        private:
            // Get the next color in the palette
            int GetNextColor()
            {
                int color = gStyle->GetColorPalette(colorIndex);
                colorIndex = (colorIndex + 1) % gStyle->GetNumberOfColors();
                return color;
            }
        };

    } // namespace plt
    */

    namespace file
    {
        class EntryCopy
        {
        private:
            typedef std::variant<
                std::unique_ptr<Int_t>,
                std::unique_ptr<Float_t>,
                std::unique_ptr<Double_t>,
                std::unique_ptr<Bool_t>,
                std::unique_ptr<Char_t>,
                std::unique_ptr<std::vector<float>>,
                std::unique_ptr<std::vector<double>>,
                std::unique_ptr<std::vector<int>>,
                std::unique_ptr<std::vector<Long64_t>>>
                supported_types;

            // Create a vector to hold pointers to the variables used for branch data and names
            std::vector<supported_types> branchData;
            std::vector<std::string> branchNames;

            static const int N = 200;
            std::vector<float> *branchData_vfloat[N];
            std::vector<double> *branchData_vdouble[N];
            std::vector<int> *branchData_vint[N];
            std::vector<Long64_t> *branchData_vint64[N];
            int vfloat_n = 0;
            int vdouble_n = 0;
            int vint_n = 0;
            int vint64_n = 0;

        public:
            EntryCopy() {}

            // Setup same branches in the destTree as the sourceTree
            void Setup(TTree *sourceTree,
                       TTree *destTree, int verbose = 0)
            {

                // Get the list of branches in the source tree
                TObjArray *branches = sourceTree->GetListOfBranches();

                // Iterate over all branches in the source tree
                for (Int_t i = 0; i < branches->GetEntries(); i++)
                {
                    TBranch *branch = (TBranch *)branches->At(i);
                    const char *branchName = branch->GetName();

                    if (sourceTree->GetBranchStatus(branchName) == false)
                        continue;

                    // Get the leaf associated with the branch
                    TLeaf *leaf = branch->GetLeaf(branchName);
                    if (!leaf)
                        continue;

                    // Determine the type of the leaf
                    const char *leafType = leaf->GetTypeName();
                    if (verbose > 0)
                        std::cout << "  Leaf type: " << leafType << " for branch " << branchName << std::endl;

                    // Create a new branch in the destination tree with the same name and type
                    if (strcmp(leafType, "Int_t") == 0)
                    {
                        auto value = std::make_unique<Int_t>();
                        destTree->Branch(branchName, value.get());
                        sourceTree->SetBranchAddress(branchName, value.get());
                        branchData.push_back(std::move(value));
                        branchNames.push_back(branchName);
                    }
                    else if (strcmp(leafType, "Float_t") == 0)
                    {
                        auto value = std::make_unique<Float_t>();
                        destTree->Branch(branchName, value.get());
                        sourceTree->SetBranchAddress(branchName, value.get());
                        branchData.push_back(std::move(value));
                        branchNames.push_back(branchName);
                    }
                    else if (strcmp(leafType, "Double_t") == 0)
                    {
                        auto value = std::make_unique<Double_t>();
                        destTree->Branch(branchName, value.get());
                        sourceTree->SetBranchAddress(branchName, value.get());
                        branchData.push_back(std::move(value));
                        branchNames.push_back(branchName);
                    }
                    else if (strcmp(leafType, "Bool_t") == 0)
                    {
                        auto value = std::make_unique<Bool_t>();
                        destTree->Branch(branchName, value.get());
                        sourceTree->SetBranchAddress(branchName, value.get());
                        branchData.push_back(std::move(value));
                        branchNames.push_back(branchName);
                    }
                    else if (strcmp(leafType, "Char_t") == 0)
                    {
                        auto value = std::make_unique<Char_t>();
                        destTree->Branch(branchName, value.get());
                        sourceTree->SetBranchAddress(branchName, value.get());
                        branchData.push_back(std::move(value));
                        branchNames.push_back(branchName);
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
                        branchNames.push_back(branchName);
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
                        branchNames.push_back(branchName);
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
                        branchNames.push_back(branchName);
                    }
                    else if (strcmp(leafType, "vector<Long64_t>") == 0)
                    {
                        auto value = new std::vector<Long64_t>();
                        branchData_vint64[vint64_n] = value;
                        destTree->Branch(branchName, &(*branchData_vint64[vint64_n]));          // Use the raw pointer
                        sourceTree->SetBranchAddress(branchName, &branchData_vint64[vint64_n]); // Use the raw pointer
                        vint64_n += 1;
                        // Convert raw pointer to unique_ptr and save to the list
                        std::unique_ptr<std::vector<Long64_t>> uniquePtr(value);
                        branchData.push_back(std::move(uniquePtr));
                        branchNames.push_back(branchName);
                    }

                    else
                    {
                        // Handle other types or arrays if necessary
                        std::cerr << "    Unsupported type: " << leafType << " for branch " << branchName << std::endl;
                        continue;
                    }
                }
            }

            // Read data from the sourceTree
            // (don't write to destTree yet,
            //  because we may want to add other things to the destTree)
            void ReadSource(TTree *sourceTree, Long64_t entry)
            {
                // Copy the data from the source tree to the destination tree
                sourceTree->GetEntry(entry);
            }

            int index(std::string key)
            {
                auto it = std::find(branchNames.begin(), branchNames.end(), key);
                int index = -1;
                if (it != branchNames.end())
                    index = std::distance(branchNames.begin(), it);
                return index;
            }

            // Get methods to easily read value with key name
            Int_t GetInt(std::string key) { return *std::get<std::unique_ptr<Int_t>>(branchData[this->index(key)]); }
            Float_t GetFloat(std::string key) { return *std::get<std::unique_ptr<Float_t>>(branchData[this->index(key)]); }
            Double_t GetDouble(std::string key) { return *std::get<std::unique_ptr<Double_t>>(branchData[this->index(key)]); }
            Bool_t GetBool(std::string key) { return *std::get<std::unique_ptr<Bool_t>>(branchData[this->index(key)]); }
            Char_t GetChar(std::string key) { return *std::get<std::unique_ptr<Char_t>>(branchData[this->index(key)]); }
            std::vector<int> *GetIntV(std::string key) { return (std::get<std::unique_ptr<std::vector<int>>>(branchData[this->index(key)])).get(); }
            std::vector<float> *GetFloatV(std::string key) { return (std::get<std::unique_ptr<std::vector<float>>>(branchData[this->index(key)])).get(); }
            std::vector<double> *GetDoubleV(std::string key) { return (std::get<std::unique_ptr<std::vector<double>>>(branchData[this->index(key)])).get(); }
        };

    } // namespace file

} // namespace iroot

#endif // util_root_hh