//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file MuAnalysis.hh
/// \brief Selection of the analysis technology

#ifndef MuAnalysis_h
#define MuAnalysis_h 1

#include <any>

#if (G4VERSION_NUMBER < 1100)
#include "g4root.hh"
#else
#include "G4AnalysisManager.hh "
#endif
#include "G4VSensitiveDetector.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"

#include "util.hh"
// #include "g4csv.hh"
// #include "g4xml.hh"

// #include "G4GenericAnalysisManager.hh"
// using G4AnalysisManager = G4GenericAnalysisManager;

namespace Analysis
{
    // Geant4 generated tuple ID
    extern int tuple_id;

    // Helper function for all analysis
    void Setup();
    bool Open(const std::string &path);
    bool Save();
    bool CreateNTuple(util::py::Dict &data, const std::string &name);
    bool FillNTuple(util::py::Dict &data);

    // // Analysis template class
    // class AnalysisTemplate
    // {
    // public:
    //     // virtual ~AnalysisTemplate();

    //     // Data dictionary
    //     util::py::Dict fdata;

    //     // Sensitive detector
    //     G4VSensitiveDetector *fDetector;
    // };

    // // Default analysis class
    // class DefaultAnalysis : public AnalysisTemplate
    // {
    // public:
    //     DefaultAnalysis();
    // };

    // Default hits/steps
    class uHit : public G4VHit
    {
    public:
        uHit();
        uHit(G4Step *step);
        // overload the allocator to use Geant4 version
        inline void *operator new(size_t);
        inline void operator delete(void *hit);

        const G4ParticleDefinition *_particle;
        int _trackID;
        int _trackPDG;
        int _parentID;
        int _parentPDG;
        double _edeposit;
        G4LorentzVector _position;
        G4LorentzVector _momentum;
    };
    extern G4Allocator<uHit> *HitAllocator;
    inline void *uHit::operator new(size_t)
    {
        return (void *)HitAllocator->MallocSingle();
    }
    inline void uHit::operator delete(void *hit)
    {
        HitAllocator->FreeSingle((uHit *)hit);
    }

    // Default hits collection
    typedef G4THitsCollection<uHit> HitsCollection;

    // Default detector
    class DefaultDetector : public G4VSensitiveDetector
    {
    public:
        DefaultDetector();
        void Initialize(G4HCofThisEvent *event) override;
        G4bool ProcessHits(G4Step *step, G4TouchableHistory *) override;
        void EndOfEvent(G4HCofThisEvent *) override;

    private:
        HitsCollection *fHitsCollection;
        util::py::Dict *fdata;
    };

}

#endif
