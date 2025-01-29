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

        static void setStyle()
        {
            gROOT->SetStyle("Modern");
            setTableauPalette();
        }
    } // namespace style

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

    namespace file
    {
        // Setup same branches in the destTree as the sourceTree
        static void CopyTreeBranches_Setup(TTree *sourceTree, TTree *destTree, const std::vector<std::string> &branchesExcluded = {})
        {
            // Get the list of branches in the source tree
            TObjArray *branches = sourceTree->GetListOfBranches();

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

                // Create a new branch in the destination tree with the same name and type
                if (strcmp(leafType, "Int_t") == 0)
                {
                    Int_t value;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "Float_t") == 0)
                {
                    Float_t value;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "Double_t") == 0)
                {
                    Double_t value;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "Bool_t") == 0)
                {
                    Bool_t value;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "Char_t") == 0)
                {
                    Char_t value;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "vector<float>") == 0)
                {
                    std::vector<float> *value = nullptr;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "vector<double>") == 0)
                {
                    std::vector<double> *value = nullptr;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "vector<int>") == 0)
                {
                    std::vector<int> *value = nullptr;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else if (strcmp(leafType, "string") == 0 || strcmp(leafType, "std::string") == 0)
                {
                    std::string *value = nullptr;
                    destTree->Branch(branchName, &value);
                    sourceTree->SetBranchAddress(branchName, &value);
                }
                else
                {
                    // Handle other types or arrays if necessary
                    std::cerr << "Unsupported type: " << leafType << " for branch " << branchName << std::endl;
                    continue;
                }
            }
        }

        // Read data from the sourceTree
        // (don't write to destTree yet,
        //  because we may want to add other things to the destTree)
        static void CopyTreeBranches_Read(TTree *sourceTree, TTree *destTree, Long64_t entry)
        {
            // Copy the data from the source tree to the destination tree
            sourceTree->GetEntry(entry);
            destTree->Fill();
        }

    } // namespace file

} // namespace iroot

#endif // util_root_hh