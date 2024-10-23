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
// #include "g4csv.hh"
// #include "g4xml.hh"

// #include "G4GenericAnalysisManager.hh"
// using G4AnalysisManager = G4GenericAnalysisManager;

namespace Analysis
{   
    // Define type of the key <single or vector>
    // Tuple can hold singles of I (int), F(float), D(double) or S(string)    
    // or vectors of  I (int), F(float), D(double)
    enum KEYTYPE
    {
        SINGLE,
        VECTOR,
    };

    // Tuple can hold data types of I (int), F(float), D(double) or S(string)    
    enum DTYPE
    {
        TUPLE_INT,
        TUPLE_FLOAT,
        TUPLE_DOUBLE,
        TUPLE_STRING,
    };


    // Hold the information of one column
    typedef struct
    {
        int key_index;
        std::string key_name;
        enum KEYTYPE key_type;
        enum KEYTYPE data_type;
        union {
            G4int data_int;
            G4float data_float;
            G4double data_doube;
            G4String *data_string;
            std::vector<G4int> *vec_int;
            std::vector<G4float> *vec_float;
            std::vector<G4double> *vec_doube;
        };
    } KEY_DATA_COLUMN_t;

    // Hold the information of everything
    typedef std::vector<KEY_DATA_COLUMN_t> KEY_DATA_t;


    void Setup();
    bool Open(const std::string &path);
    bool Save();
    bool CreateNTuple(KEY_DATA_t &data, const std::string &name);
    bool FillNTuple(KEY_DATA_t &data);

}

#endif
