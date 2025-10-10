//________________________________________________________________________________________
/*!

\program gevgen_flux

\brief   A GENIE atmospheric neutrino event generation application.

         *** Synopsis :

           gevgen_flux [-h]
                       [-r run#]
                        -f flux
                        -g geometry
                       [-R coordinate_rotation_matrix]
                       [-t geometry_top_volume_name]
                       [-m max_path_lengths_xml_file]
                       [-L geometry_length_units]
                       [-D geometry_density_units]
                       <-n n_of_events,
                        -e exposure_in_kton_x_yrs
                        -T exposure_in_seconds >
                        -E min_energy,max_energy
                       [-o output_event_file_prefix]
                       [--flux-ray-generation-surface-distance ]
                       [--flux-ray-generation-surface-radius   ]
                       [--seed random_number_seed]
                       [--cross-sections xml_file]
                       [--event-generator-list list_name]
                       [--tune genie_tune]
                       [--message-thresholds xml_file]
                       [--unphysical-event-mask mask]
                       [--event-record-print-level level]
                       [--mc-job-status-refresh-rate  rate]
                       [--cache-file root_file]

         *** Options :

           [] Denotes an optional argument.
           <> Denotes a set of arguments out of which only one can be set.

           -h
              Prints out the syntax and exits
           -r
              Specifies the MC run number
              [default: 100000000]
           -f
              Specifies the input flux files
              A set of histograms stored in a ROOT file.
                  The general syntax is:
                      -f /path/histogram_file.root,neutrino_code[histo_name],...
                  [Notes]
                  - The neutrino codes are the PDG ones.
                  - The 'neutrino_code[histogram_name]' part of the option can be
                    repeated multiple times (separated by commas), once for each
                    flux neutrino species you want to consider, eg
                    '-f somefile.root,12[nuehst],-12[nuebarhst],14[numuhst]'
                  - Code implicitly assumes the binning for multiple flavors
                    is the same for all histograms (no checks are made)
                  - Note that the relative normalization of the flux histograms
                    is taken into account and is reflected in the relative
                    frequency of flux neutrinos thrown by the flux driver
                  - Variable bin width flux histograms are modified to account
                    for incorrect sampling of TH1D::GetRandom() prior to ROOT
                    version 9.99.99
                    + notation such as "12[nuehst]WIDTH" forces bin width
                      adjustment no matter variable bin width histogram or not
                    + notation such as "12[nuehst]NOWIDTH" prevents bin width
                      adjustment
                  - When using flux from histograms then there is no point in using
                    a 'detailed detector geometry description' as your flux input
                    contains no directional information for those flux neutrinos.
                    The neutrino direction is conventionally set to be +z {x=0,y=0}.
                    So, when using this option you must be using a simple 'target mix'
                    See the -g option for possible geometry settings.
                    If you want to use the detailed detector geometry description
                    then you should be feeding this driver with the JNUBEAM flux
                    simulation outputs.
                  - When using flux from histograms there is no branch with neutrino
                    parent information (code, decay mode, 4-momentum at prod & decay)
                    added in the output event tree as your flux input contains no
                    such information. If you want to be getting the neutrino parent
                    information written out as well then you should be feeding this
                    driver with the JNUBEAM flux simulation outputs.
                  - Note that the relative normalization of the flux histograms is
                    taken into account and is reflected in the relative frequency
                    of flux neutrinos thrown by the flux driver
                  [Examples]
                  - To use the histogram 'h100' (representing the nu_mu flux) and
                    the histogram 'h300' (representing the nu_e flux) and the
                    histogram 'h301' (representing the nu_e_bar flux) from the
                    flux.root file in /path/
                    type:
                      '-f /path/flux.root,14[h100],12[h300],-12[h301]                      
           -g
              Input 'geometry'.
              This option can be used to specify any of:
              1 > A ROOT file containing a ROOT/GEANT geometry description
                  [Examples]
                  - To use the master volume from the ROOT geometry stored
                    in the nd280-geom.root file, type:
                    '-g /some/path/nd280-geom.root'
              2 > A mix of target materials, each with its corresponding weight,
                  typed as a comma-separated list of nuclear pdg codes (in the
                  std PDG2006 convention: 10LZZZAAAI) with the weight fractions
                  in brackets, eg code1[fraction1],code2[fraction2],...
                  If that option is used (no detailed input geometry description)
                  then the interaction vertices are distributed in the detector
                  by the detector MC.
                  [Examples]
                  - To use a target mix of 89% O16 and 11% H, type:
                    '-g 1000080160[0.89],1000010010[0.11]'
                  - To use a target which is 100% C12, type:
                    '-g 1000060120'
           -R
              Input rotation matrix for transforming the flux neutrino coordinates
              from the default Topocentric Horizontal (see GENIE manual) coordinate
              system to the user-defined topocentric coordinate system.
              The rotation is specified by the 3 Euler angles (phi, theta, psi).
              The Euler angles are input as a comma separated list as:
              `-R <convention>:phi,theta,psi',
              where <convention> is either X (for X-convention), Y (for Y-convention),
              X^-1 or Y^-1 (as previously, but using the inverse matrix).
              By default, the X-convention (rotation about Z-axis, then about the
              new X-axis, then about the Z-axis) is used.
              Notes:
              - (Extract from TRotation documentation)
               "Euler angles usually define the rotation of the new coordinate
                system with respect to the original system, however, the TRotation
                class specifies the rotation of the object in the original system
                (an active rotation). To recover the usual Euler rotations (ie.
                rotate the system not the object), you must take the inverse of
                the rotation."
              Examples:
              1. To set the Euler angles phi=3.14, theta=1.28, psi=1.0 using the
                 X-convention, type: `-R 3.14,1.28,1.0', or `-R X:3.14,1.28,1.0'
              2. To set the Euler angles phi=3.14, theta=1.28, psi=1.0 using the
                 Y-convention, type: `-R Y:3.14,1.28,1.0'
              3. To set the Euler angles phi=3.14, theta=1.28, psi=1.0 using the
                 Y-convention, and then use the inverse rotation matrix, type:
                 `-R Y^-1:3.14,1.28,1.0'
           -L
              Input geometry length units, eg 'm', 'cm', 'mm', ...
              [default: 'mm']
           -D
              Input geometry density units, eg 'g_cm3', 'clhep_def_density_unit',...
              [default: 'g_cm3']
           -t
              Input 'top volume' for event generation -
              can be used to force event generation in given sub-detector
              [default: the 'master volume' of the input geometry]
              You can also use the -t option to switch generation on/off at
              multiple volumes as, for example, in:
              `-t +Vol1-Vol2+Vol3-Vol4',
              `-t "+Vol1 -Vol2 +Vol3 -Vol4"',
              `-t -Vol2-Vol4+Vol1+Vol3',
              `-t "-Vol2 -Vol4 +Vol1 +Vol3"'
              where:
              "+Vol1" and "+Vol3" tells GENIE to `switch on'  Vol1 and Vol3, while
              "-Vol2" and "-Vol4" tells GENIE to `switch off' Vol2 and Vol4.
              If the very first character is a '+', GENIE will neglect all volumes
              except the ones explicitly turned on. Vice versa, if the very first
              character is a `-', GENIE will keep all volumes except the ones
              explicitly turned off (feature contributed by J.Holeczek).
           -n
              Specifies how many events to generate.
           -e
              Specifies requested exposure in terms of kton*yrs.
           -T
              Specifies requested exposure in terms of seconds.
           -E
              Specifies the neutrino energy in GeV.
              Must be a comma-separated pair of numbers, eg `-E 0.3,70'
              [default: 0.5,50]
           --flux-ray-generation-surface-distance
           --flux-ray-generation-surface-radius
              See the User & Physics Manual for a graphical representation of the flux
              ray generation surface: For a given zenith \theta and azimuthal angle \phi,
              flux rays are produced within the area of a circle of radius Rt,
              which is tangetial to a point P on a sphere of radius Rl, centred at the
              detector. The point P has polar coordinates \theta and \phi.
              The argument --flux-ray-generation-surface-distance sets Rl, while
              the argument --flux-ray-generation-surface-distance sets Rt.
              SI units are used.
           -o
              Sets the prefix of the output event file.
              The output filename is built as:
              [prefix].[run_number].[event_tree_format].[file_format]
              The default output filename is:
              gntp.[run_number].ghep.root
              This cmd line arguments lets you override 'gntp'
           --seed
              Random number seed.
           --cross-sections
              Name (incl. full path) of an XML file with pre-computed
              cross-section values used for constructing splines.
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
              [default: "Default"].
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
              Multiple files, delimited with a `:' can be specified.
           --unphysical-event-mask
              Allows users to specify a 16-bit mask to allow certain types of
              unphysical events to be written in the output file.
              [default: all unphysical events are rejected]
           --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().
           --mc-job-status-refresh-rate
              Allows users to customize the refresh rate of the status file.
           --cache-file
              Allows users to specify a cache file so that the cache can be
              re-used in subsequent MC jobs.

         *** Examples:

           (1) Generate 100k events (run number 999210) in the energy range 1-10 GeV
               for nu_e and nu_mu only, using the sdave_numu07.dat FLUKA flux file for
               nu_mu and the sdave_nue07.dat file for nu_e (files in /data/flx/).
               Use the detector geometry in the /data/geo/SuperK.root file, where the
               geometry length and density units are m and kgr/m^3. Generate events over
               the entire geometry volume. Pre-computed cross-section data are loaded
               from /data/xsec.xml.

               % gevgen_flux -r 999210 -n 100000 -E 1,10
                       -f FLUKA:/data/flx/sdave_numu07.dat[14],/data/flx/sdave_nue07.dat[12]
                       -g /data/geo/SuperK.root -L "m" -D "kg_m3"
                       --cross-sections /data/xsec.xml

           (2) Like above but, instead of generating events in a realistic detector
               geometry, use a simple target mix (88.79% O16 + 11.21% H, i.e. `water')

               % gevgen_flux -r 999210 -n 100000 -E 1,10
                       -f /data/flux/sdave_numu07.dat[14],/data/flux/sdave_nue07.dat[12]
                       -g 1000080160[0.8879],1000010010[0.1121]
                       --cross-sections /data/xsec.xml

                ... to add more

         Please read the GENIE User Manual for more information.

\created August 20, 2010

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
         University of Liverpool

         Torben Ferber <torben.ferber \at DESY.DE>
         DESY

         Hugh Gallagher <hugh.gallagher \at tufts.edu>
         Tufts University

         Tarak Thakore <tarak \at mailhost.tifr.res.in>
         Tata Institute of Fundamental Research

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <iomanip>
#include <cmath>

#include <TRotation.h>
#include <TMath.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#include <TSystem.h>

#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/EventGen/GMCJDriver.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "Tools/Flux/GFLUKAAtmoFlux.h"
#include "Tools/Flux/GBGLRSAtmoFlux.h"
#include "Tools/Flux/GHAKKMAtmoFlux.h"
#include "Tools/Flux/GCylindTH1Flux.h"
#endif

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#include "Tools/Geometry/GeoUtils.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"
#include "Tools/Geometry/PointGeomAnalyzer.h"
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;
using std::setprecision;

using namespace genie;
using namespace genie::flux;

void            GetCommandLineArgs (int argc, char ** argv);
void            PrintSyntax        (void);
GFluxI*      GetFluxHistDriver(void);
GeomAnalyzerI * GetGeometry        (void);
float total_flux_input = 0;

// User-specified options:
//
Long_t          gOptRunNu;                     // run number
string          gOptFluxSim;                   // flux simulation (FLUKA, BGLRS or HAKKM)
string          gOptFluxFile;                  // ROOT file with JNUBEAM flux ntuple
map<int,TH1D*>  gOptFluxHst;                   // flux histos (nu pdg  -> spectrum)  / if not using beam sim ntuples
bool            gOptUsingRootGeom = false;     // using root geom or target mix?
map<int,double> gOptTgtMix;                    // target mix  (tgt pdg -> wght frac) / if not using detailed root geom
string          gOptRootGeom;                  // input ROOT file with realistic detector geometry
string          gOptRootGeomTopVol = "";       // input geometry top event generation volume
double          gOptGeomLUnits = 0;            // input geometry length units
double          gOptGeomDUnits = 0;            // input geometry density units
string          gOptExtMaxPlXml;               // max path lengths XML file for input geometry
int             gOptNev = -1;                  // exposure - in terms of number of events
double          gOptKtonYrExposure = -1;       // exposure - in terms of kton*yrs
double          gOptSecExposure = -1;          // exposure - in terms of seconds
double          gOptEvMin;                     // minimum neutrino energy
double          gOptEvMax;                     // maximum neutrino energy
string          gOptEvFilePrefix;              // event file prefix
TRotation       gOptRot;                       // coordinate rotation matrix: topocentric horizontal -> user-defined topocentric system
long int        gOptRanSeed;                   // random number seed
string          gOptInpXSecFile;               // cross-section splines
double          gOptRL = -1;                   // distance of flux ray generation surface (m)
double          gOptRT = -1;                   // radius of flux ray generation surface (m)

// Defaults:
//
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // def event tree format
string          kDefOptEvFilePrefix = "gntp";  // def output prefix (override with -o)
string          kDefOptGeomLUnits   = "mm";    // def geom length units (override with -L)
string          kDefOptGeomDUnits   = "g_cm3"; // def geom density units (override with -D)
double          kDefOptEvMin        =  0.5;    // min neutrino energy (override with -E)
double          kDefOptEvMax        = 50.0;    // max neutrino energy (override with -E)

//________________________________________________________________________________________
int main(int argc, char** argv)
{
  GFluxI* atmo_flux_driver;
  double total_flux, expected_neutrinos;

  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gmkspl", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // Iinitialization of random number generators, cross-section table, messenger, cache etc...
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, true);

  // get geometry driver
  GeomAnalyzerI * geom_driver = GetGeometry();

  if (gOptRT < 0) {
    gOptRT = 1000; // m
    LOG("gevgen_flux", pWARN) << "Warning! Flux surface transverse radius not specified so using default value of " << gOptRT << " meters!";
  }

  if (gOptRL < 0) {
    gOptRL = 1000; // m
    LOG("gevgen_flux", pWARN) << "Warning! Flux surface longitudinal radius not specified so using default value of " << gOptRL << " meters!";
  }

  // get flux driver
  atmo_flux_driver = GetFluxHistDriver();

  // Cast to GFluxI, the generic flux driver interface
  GFluxI *flux_driver = dynamic_cast<GFluxI *>(atmo_flux_driver);

  // create the GENIE monte carlo job driver
  GMCJDriver* mcj_driver = new GMCJDriver;
  mcj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();
  /* Note: For the method of calculating the total number of events using a
   * livetime we *need* to force a single probability scale. So if you change
   * the next line you need to make sure that the user didn't specify the -T
   * option. */
  mcj_driver->ForceSingleProbScale();

  // initialize an ntuple writer
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // Create a MC job monitor for a periodically updated status file
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  total_flux = total_flux_input; //atmo_flux_driver->GetTotalFluxInEnergyRange();
  LOG("gevgen_flux", pNOTICE) << "Total neutrino flux is " << setprecision(2) << total_flux << " neutrinos per second.";
  if (gOptSecExposure > 0) {
    /* Calculate the expected value of the total number of neutrinos we need to
     * throw. We do this by multiplying the total flux by the exposure time in
     * seconds and the area of the flux surface. */
    expected_neutrinos = total_flux*gOptSecExposure; // *atmo_flux_driver->GetFluxSurfaceArea();
    LOG("gevgen_flux", pNOTICE) << "Simulating an exposure of " << setprecision(0) << gOptSecExposure << " seconds which corresponds to a total of " << setprecision(0) << expected_neutrinos << " neutrinos.";
  }

  // event loop
  for(int iev = 0; gOptNev > 0 ? iev < gOptNev : 1; iev++) {

    // generate next event
    EventRecord* event = mcj_driver->GenerateEvent();

    // print-out
    LOG("gevgen_flux", pNOTICE) << "Generated event: " << *event;

    // save the event, refresh the mc job monitor
    ntpw.AddEventRecord(iev, event);
    mcjmonitor.Update(iev,event);

    // clean-up
    delete event;

    if (gOptSecExposure > 0 && mcj_driver->NFluxNeutrinos()/mcj_driver->GlobProbScale() > expected_neutrinos) {
      break;
    }
  }

  // save the event file
  ntpw.Save();

  // clean-up
  delete geom_driver;
  delete atmo_flux_driver;
  delete mcj_driver;

  return 0;
}
//________________________________________________________________________________________
GeomAnalyzerI* GetGeometry(void)
{
  GeomAnalyzerI * geom_driver = 0;

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__

  if(gOptUsingRootGeom) {
    //
    // *** Using a realistic root-based detector geometry description
    //

    // creating & configuring a root geometry driver
    geometry::ROOTGeomAnalyzer * rgeom =
            new geometry::ROOTGeomAnalyzer(gOptRootGeom);
    rgeom -> SetLengthUnits  (gOptGeomLUnits);
    rgeom -> SetDensityUnits (gOptGeomDUnits);
    rgeom -> SetTopVolName   (gOptRootGeomTopVol);
    // getting the bounding box dimensions along z so as to set the
    // appropriate upstream generation surface for the JPARC flux driver
    TGeoVolume * topvol = rgeom->GetGeometry()->GetTopVolume();
    if(!topvol) {
      LOG("gevgen_flux", pFATAL) << " ** Null top ROOT geometry volume!";
      gAbortingInErr = true;
      exit(1);
    }

    /* If flux generation surface isn't defined, get the bounding box for the
     * geometry and set something appropriate. */
    TGeoShape *bounding_box = topvol->GetShape();
    TGeoBBox *box = (TGeoBBox *) bounding_box;
    double dx = box->GetDX()*rgeom->LengthUnits();
    double dy = box->GetDY()*rgeom->LengthUnits();
    double dz = box->GetDZ()*rgeom->LengthUnits();

    if (gOptRL < 0 && gOptRT < 0) {
      gOptRL = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
      gOptRT = gOptRL;
      LOG("gevgen_flux", pNOTICE) << "Setting flux longitudinal and transverse radius to " << setprecision(2) << gOptRL << " meters based on bounding box of ROOT geometry.";
    }

    // switch on/off volumes as requested
    if ( (gOptRootGeomTopVol[0] == '+') || (gOptRootGeomTopVol[0] == '-') ) {
      bool exhaust = (*gOptRootGeomTopVol.c_str() == '+');
      utils::geometry::RecursiveExhaust(topvol, gOptRootGeomTopVol, exhaust);
    }

    // casting to the GENIE geometry driver interface
    geom_driver = dynamic_cast<GeomAnalyzerI *> (rgeom);
  }
  else {
    //
    // *** Using a 'point' geometry with the specified target mix
    // *** ( = a list of targets with their corresponding weight fraction)
    //

    // creating & configuring a point geometry driver
    geometry::PointGeomAnalyzer * pgeom =
              new geometry::PointGeomAnalyzer(gOptTgtMix);
    // casting to the GENIE geometry driver interface
    geom_driver = dynamic_cast<GeomAnalyzerI *> (pgeom);
  }

#else
  LOG("gevgen_flux", pFATAL) << "You need to enable the GENIE geometry drivers first!";
  LOG("gevgen_flux", pFATAL) << "Use --enable-geom-drivers at the configuration step.";
  gAbortingInErr = true;
  exit(1);
#endif

  return geom_driver;
}
//________________________________________________________________________________________
GFluxI* GetFluxHistDriver(void)
{
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__

    //
    // *** Using fluxes from histograms (for all specified neutrino species)
    //
    GFluxI * flux_driver = 0;
    flux::GCylindTH1Flux * hst_flux_driver   = 0;

    // creating & configuring a generic GCylindTH1Flux flux driver
    TVector3 bdir (-85, 0, -(85.47+5)); // direction
    TVector3 bspot(90/5, 0, 85.47/5); // beam spot, After trial and error, it turns out that the unit is [m]

    // TVector3 bdir (0,0,-1); // dir along +z
    // TVector3 bspot(0,0,85.47); // beam spot, After trial and error, it turns out that the unit is [m]    
    hst_flux_driver = new flux::GCylindTH1Flux;
    hst_flux_driver->SetNuDirection      (bdir);
    hst_flux_driver->SetBeamSpot         (bspot);
    hst_flux_driver->SetTransverseRadius (1);
    map<int,TH1D*>::iterator it = gOptFluxHst.begin();
    for( ; it != gOptFluxHst.end(); ++it) {
        int    pdg_code = it->first;
        TH1D * spectrum = new TH1D(*(it->second));
        hst_flux_driver->AddEnergySpectrum(pdg_code, spectrum);

        total_flux_input += spectrum->Integral("width"); // Total flux is sum of integrals of all spectra
    }
    // casting to the GENIE flux driver interface
    flux_driver = dynamic_cast<GFluxI *> (hst_flux_driver);

#else
  LOG("gevgen_flux", pFATAL) << "You need to enable the GENIE flux drivers first!";
  LOG("gevgen_flux", pFATAL) << "Use --enable-flux-drivers at the configuration step.";
  gAbortingInErr = true;
  exit(1);
#endif

  return flux_driver;
}
//________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
// Get the command line arguments

  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  LOG("gevgen_flux", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  //
  // *** run number
  //
  if( parser.OptionExists('r') ) {
    LOG("gevgen_flux", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_flux", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 100000000;
  } //-r

  //
  // *** exposure
  //

  // in number of events
  bool have_required_statistics = false;
  if( parser.OptionExists('n') ) {
    LOG("gevgen_flux", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
    have_required_statistics = true;
  }//-n?
  // or, in kton*yrs
  if( parser.OptionExists('e') ) {
    if(have_required_statistics) {
      LOG("gevgen_flux", pFATAL)
         << "Can't request exposure both in terms of number of events and  kton*yrs"
         << "\nUse just one of the -n and -e options";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
    LOG("gevgen_flux", pDEBUG)
        << "Reading requested exposure in kton*yrs";
    gOptKtonYrExposure = parser.ArgAsDouble('e');
    have_required_statistics = true;
  }//-e?

  if (parser.OptionExists('T')) {
    if (have_required_statistics) {
      LOG("gevgen_flux", pFATAL)
         << "Can't request exposure both in terms of number of events or kton*yrs and time"
         << "\nUse just one of the -n, -e, or -T options";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
    LOG("gevgen_flux", pDEBUG)
        << "Reading requested exposure in seconds";
    gOptSecExposure = parser.ArgAsDouble('T');
    have_required_statistics = true;
  }

  if(!have_required_statistics) {
    LOG("gevgen_flux", pFATAL)
       << "You must request exposure either in terms of number of events and  kton*yrs"
       << "\nUse any of the -n, -e options";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  //
  // *** event file prefix
  //
  if( parser.OptionExists('o') ) {
    LOG("gevgen_flux", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_flux", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  //
  // *** neutrino energy range
  //
  if( parser.OptionExists('E') ) {
    LOG("gevgen_flux", pINFO) << "Reading neutrino energy range";
    string nue = parser.ArgAsString('E');

    // must be a comma separated set of values
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0.);
       gOptEvMin = emin;
       gOptEvMax = emax;
    } else {
      LOG("gevgen_flux", pFATAL)
        << "Invalid energy range. Use `-E emin,emax', eg `-E 0.5,100.";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
  } else {
     LOG("gevgen_flux", pNOTICE)
        << "No -e option. Using default energy range";
     gOptEvMin = kDefOptEvMax;
     gOptEvMax = kDefOptEvMax;
  }

  //
  // *** flux files
  //
  // syntax:
  // simulation:/path/file.data[neutrino_code],/path/file.data[neutrino_code],...
  //
  if( parser.OptionExists('f') ) {
    LOG("gevgen_flux", pDEBUG) << "Getting input flux files";
    string flux = parser.ArgAsString('f');


   // Using flux from histograms
   // Extract the root file name & the list of histogram names & neutrino
   // species (specified as 'filename,histo1[species1],histo2[species2],...')
   // for variable width histograms, default to multiply in the width
   //   histo1[species1]WIDTH = multiply in the width
   //   histo1[species1]NOWIDTH = don't multiply in the width
   // See documentation on top section of this file.
   //
   vector<string> fluxv = utils::str::Split(flux,",");
   if(fluxv.size()<2) {
      LOG("gevgen_flux", pFATAL)
         << "You need to specify both a flux ntuple ROOT file "
         << " _AND_ a detector location";
      PrintSyntax();
      exit(1);
   }
   gOptFluxFile = fluxv[0];
   bool accessible_flux_file =
         !(gSystem->AccessPathName(gOptFluxFile.c_str()));
   if (!accessible_flux_file) {
      LOG("gevgen_flux", pFATAL)
         << "Can not access flux file: " << gOptFluxFile;
      PrintSyntax();
      exit(1);
   }
   // Extract energy spectra for all specified neutrino species
   TFile flux_file(gOptFluxFile.c_str(), "read");
   for(unsigned int inu=1; inu<fluxv.size(); inu++) {
      string nutype_and_histo = fluxv[inu];
      string::size_type open_bracket  = nutype_and_histo.find("[");
      string::size_type close_bracket = nutype_and_histo.find("]");
      if (open_bracket ==string::npos ||
            close_bracket==string::npos)
      {
            LOG("gevgen_flux", pFATAL)
               << "You made an error in specifying the flux histograms";
            PrintSyntax();
            exit(1);
      }
      string::size_type ibeg = 0;
      string::size_type iend = open_bracket;
      string::size_type jbeg = open_bracket+1;
      string::size_type jend = close_bracket;
      string nutype = nutype_and_histo.substr(ibeg,iend-ibeg);
      string histo  = nutype_and_histo.substr(jbeg,jend-jbeg);
      string extra  = nutype_and_histo.substr(jend+1,string::npos);
      std::transform(extra.begin(),extra.end(),extra.begin(),::toupper);
      // access specified histogram from the input root file
      TH1D * ihst = (TH1D*) flux_file.Get(histo.c_str());
      if(!ihst) {
            LOG("gevgen_flux", pFATAL)
            << "Can not find histogram: " << histo
            << " in flux file: " << gOptFluxFile;
            PrintSyntax();
            exit(1);
      }
      // create a local copy of the input histogram
      TString origname = ihst->GetName();
      TString tmpname; tmpname.Form("%s_", origname.Data());
      // Copy in the flux histogram from the root file
      // use Clone rather than assuming fix bin widths and rebooking
      TH1D* spectrum = (TH1D*)ihst->Clone();
      spectrum->SetNameTitle(tmpname.Data(),ihst->GetName());
      spectrum->SetDirectory(0);

      // get rid of original
      delete ihst;
      // rename copy
      spectrum->SetName(origname.Data());

      bool force_binwidth = false;
#if ROOT_VERSION_CODE <= ROOT_VERSION(9,99,99)
      // GetRandom() sampling on variable bin width histograms does not
      // correctly account for bin widths for all versions of ROOT prior
      // to (currently forever).  At some point this might change and
      // the necessity of this code snippet will go away
      TAxis* xaxis = spectrum->GetXaxis();
      if (xaxis->IsVariableBinSize()) force_binwidth = true;
#endif
      if ( extra == "WIDTH"   ) force_binwidth = true;
      if ( extra == "NOWIDTH" ) force_binwidth = false;
      if ( force_binwidth ) {
         LOG("gevgen_flux", pNOTICE)
            << "multiplying by bin width for histogram " << histo
            << " as " << spectrum->GetName() << " for nutype " << nutype
            << " from " << gOptFluxFile;
         for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ++ibin) {
            double data = spectrum->GetBinContent(ibin);
            double width = spectrum->GetBinWidth(ibin);
            spectrum->SetBinContent(ibin,data*width);
         }
      }

      // convert neutrino name -> pdg code
      int pdg = atoi(nutype.c_str());
      if(!pdg::IsNeutrino(pdg) && !pdg::IsAntiNeutrino(pdg)) {
            LOG("gevgen_flux", pFATAL)
               << "Unknown neutrino type: " << nutype;
            PrintSyntax();
            exit(1);
      }
      // store flux neutrino code / energy spectrum
      LOG("gevgen_flux", pDEBUG)
         << "Adding energy spectrum for flux neutrino: pdg = " << pdg;
      gOptFluxHst.insert(map<int, TH1D*>::value_type(pdg, spectrum));
   }//inu
   if(gOptFluxHst.size()<1) {
      LOG("gevgen_flux", pFATAL)
         << "You have not specified any flux histogram!";
      PrintSyntax();
      exit(1);
   }
   flux_file.Close();
   } else {
    LOG("gevgen_flux", pFATAL) << "No flux info was specified! Use the -f option.";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // *** options to fine tune the flux ray generation surface

  if( parser.OptionExists("flux-ray-generation-surface-distance") ) {
    LOG("gevgen_flux", pINFO)
      << "Reading distance of flux ray generation surface";
    gOptRL = parser.ArgAsDouble("flux-ray-generation-surface-distance");
  } else {
    LOG("gevgen_flux", pINFO)
      << "Unspecified distance of flux ray generation surface - Using default";
  }

  if( parser.OptionExists("flux-ray-generation-surface-radius") ) {
    LOG("gevgen_flux", pINFO)
      << "Reading radius of flux ray generation surface";
    gOptRT = parser.ArgAsDouble("flux-ray-generation-surface-radius");
  } else {
    LOG("gevgen_flux", pINFO)
      << "Unspecified radius of flux ray generation surface - Using default";
  }

  //
  // *** geometry
  //
  string geom = "";
  string lunits, dunits;
  if( parser.OptionExists('g') ) {
    LOG("gevgen_flux", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
        utils::system::FileExists(geom.c_str());
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    }
  } else {
      LOG("gevgen_flux", pFATAL)
        << "No geometry option specified - Exiting";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // legth units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_flux", pDEBUG)
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_flux", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // density units:
     if( parser.OptionExists('D') ) {
        LOG("gevgen_flux", pDEBUG)
           << "Checking for input geometry density units";
        dunits = parser.ArgAsString('D');
     } else {
        LOG("gevgen_flux", pDEBUG) << "Using default geometry density units";
        dunits = kDefOptGeomDUnits;
     } // -D
     gOptGeomLUnits = genie::utils::units::UnitFromString(lunits);
     gOptGeomDUnits = genie::utils::units::UnitFromString(dunits);

     // check whether an event generation volume name has been
     // specified -- default is the 'top volume'
     if( parser.OptionExists('t') ) {
        LOG("gevgen_flux", pDEBUG) << "Checking for input volume name";
        gOptRootGeomTopVol = parser.ArgAsString('t');
     } else {
        LOG("gevgen_flux", pDEBUG) << "Using the <master volume>";
     } // -t

     // check whether an XML file with the maximum (density weighted)
     // path lengths for each detector material is specified -
     // otherwise will compute the max path lengths at job init
     if( parser.OptionExists('m') ) {
        LOG("gevgen_flux", pDEBUG)
              << "Checking for maximum path lengths XML file";
        gOptExtMaxPlXml = parser.ArgAsString('m');
     } else {
        LOG("gevgen_flux", pDEBUG)
           << "Will compute the maximum path lengths at job init";
        gOptExtMaxPlXml = "";
     } // -m
  } // using root geom?

  else {
    // User has specified a target mix.
    // Decode the list of target pdf codes & their corresponding weight fraction
    // (specified as 'pdg_code_1[fraction_1],pdg_code_2[fraction_2],...')
    // See documentation on top section of this file.
    //
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(geom,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
    } else {
      vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
      for( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
         string tgt_with_wgt = *tgtmix_iter;
         string::size_type open_bracket  = tgt_with_wgt.find("[");
         string::size_type close_bracket = tgt_with_wgt.find("]");
         if (open_bracket ==string::npos ||
             close_bracket==string::npos)
         {
             LOG("gevgen_flux", pFATAL)
                << "You made an error in specifying the target mix";
             PrintSyntax();
             gAbortingInErr = true;
             exit(1);
         }
         string::size_type ibeg = 0;
         string::size_type iend = open_bracket;
         string::size_type jbeg = open_bracket+1;
         string::size_type jend = close_bracket;
         int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend-ibeg).c_str());
         double wgt = atof(tgt_with_wgt.substr(jbeg,jend-jbeg).c_str());
         LOG("gevgen_flux", pDEBUG)
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));

      }// tgtmix_iter
    } // >1 materials in mix
  } // using tgt mix?

  //
  // Coordinate rotation matrix
  //
  gOptRot.SetToIdentity();
  if( parser.OptionExists('R') ) {
    string rotarg = parser.ArgAsString('R');
    //get convention
    string::size_type j = rotarg.find_first_of(":",0);
    string convention = "";
    if(j==string::npos) { convention = "X"; }
    else                { convention = rotarg.substr(0,j); }
    //get angles phi,theta,psi
    rotarg.erase(0,j+1);
    vector<string> euler_angles = utils::str::Split(rotarg,",");
    if(euler_angles.size() != 3) {
       LOG("gevgen_flux", pFATAL)
         << "You didn't specify all 3 Euler angles using the -R option";
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }
    double phi   = atof(euler_angles[0].c_str());
    double theta = atof(euler_angles[1].c_str());
    double psi   = atof(euler_angles[2].c_str());
    //set Euler angles using appropriate convention
    if(convention.find("X")!=string::npos ||
       convention.find("x")!=string::npos)
    {
       LOG("gevgen_flux", pNOTICE) << "Using X-convention for input Euler angles";
       gOptRot.SetXEulerAngles(phi,theta,psi);
    } else
    if(convention.find("Y")!=string::npos ||
       convention.find("y")!=string::npos)
    {
       LOG("gevgen_flux", pNOTICE) << "Using Y-convention for input Euler angles";
       gOptRot.SetYEulerAngles(phi,theta,psi);
    } else {
       LOG("gevgen_flux", pFATAL)
         << "Unknown Euler angle convention. Please use the X- or Y-convention";
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }
    //invert?
    if(convention.find("^-1")!=string::npos) {
       LOG("gevgen_flux", pNOTICE) << "Inverting rotation matrix";
       gOptRot.Invert();
    }
  }

  //
  // *** random number seed
  //
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_flux", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_flux", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  //
  // *** input cross-section file
  //
  if( parser.OptionExists("cross-sections") ) {
    LOG("gevgen_flux", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("cross-sections");
  } else {
    LOG("gevgen_flux", pINFO) << "Unspecified cross-section file";
    gOptInpXSecFile = "";
  }

  //
  // print-out summary
  //

  PDGLibrary * pdglib = PDGLibrary::Instance();

  ostringstream gminfo;
  if (gOptUsingRootGeom) {
    gminfo << "Using ROOT geometry - file: " << gOptRootGeom
           << ", top volume: "
           << ((gOptRootGeomTopVol.size()==0) ? "<master volume>" : gOptRootGeomTopVol)
           << ", max{PL} file: "
           << ((gOptExtMaxPlXml.size()==0) ? "<none>" : gOptExtMaxPlXml)
           << ", length  units: " << lunits
           << ", density units: " << dunits;
  } else {
    gminfo << "Using target mix - ";
    map<int,double>::const_iterator iter;
    for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
          int    pdg_code = iter->first;
          double wgt      = iter->second;
          TParticlePDG * p = pdglib->Find(pdg_code);
          if(p) {
            string name = p->GetName();
            gminfo << "(" << name << ") -> " << 100*wgt << "% / ";
          }//p?
    }
  }

  ostringstream fluxinfo;
  fluxinfo << "Using " << gOptFluxSim << " flux files: ";


  ostringstream expinfo;
  if(gOptNev > 0)            { expinfo << gOptNev            << " events";   }
  if(gOptKtonYrExposure > 0) { expinfo << gOptKtonYrExposure << " kton*yrs"; }

  ostringstream rotation;
  rotation << "\t| " <<  gOptRot.XX() << "  " << gOptRot.XY() << "  " << gOptRot.XZ() << " |\n";
  rotation << "\t| " <<  gOptRot.YX() << "  " << gOptRot.YY() << "  " << gOptRot.YZ() << " |\n";
  rotation << "\t| " <<  gOptRot.ZX() << "  " << gOptRot.ZY() << "  " << gOptRot.ZZ() << " |\n";

  LOG("gevgen_flux", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_flux job configuration");

  LOG("gevgen_flux", pNOTICE)
     << "\n"
     << "\n @@ Run number: " << gOptRunNu
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Using cross-section file: " << gOptInpXSecFile
     << "\n @@ Geometry"
     << "\n\t" << gminfo.str()
     << "\n @@ Flux"
     << "\n\t" << fluxinfo.str()
     << "\n @@ Exposure"
     << "\n\t" << expinfo.str()
     << "\n @@ Cuts"
     << "\n\t Using energy range = (" << gOptEvMin << " GeV, " << gOptEvMax << " GeV)"
     << "\n @@ Coordinate transformation (Rotation THZ -> User-defined coordinate system)"
     << "\n" << rotation.str()
     << "\n\n";

  //
  // final checks
  //
  if(gOptKtonYrExposure > 0) {
    LOG("gevgen_flux", pFATAL)
      << "\n Option to set exposure in terms of kton*yrs not supported just yet!"
      << "\n Try the -n option instead";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }
}
//________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_flux", pFATAL)
   << "\n **Syntax**"
   << "\n gevgen_flux [-h]"
   << "\n           [-r run#]"
   << "\n            -f simulation:flux_file[neutrino_code],..."
   << "\n            -g geometry"
   << "\n           [-R coordinate_rotation_matrix]"
   << "\n           [-t geometry_top_volume_name]"
   << "\n           [-m max_path_lengths_xml_file]"
   << "\n           [-L geometry_length_units]"
   << "\n           [-D geometry_density_units]"
   << "\n           <-n n_of_events,"
   << "\n            -e exposure_in_kton_x_yrs"
   << "\n            -T exposure_in_seconds>"
   << "\n            -E min_energy,max_energy"
   << "\n           [-o output_event_file_prefix]"
   << "\n           [--flux-ray-generation-surface-distance]"
   << "\n           [--flux-ray-generation-surface-radius]"
   << "\n           [--seed random_number_seed]"
   << "\n            --cross-sections xml_file"
   << "\n           [--event-generator-list list_name]"
   << "\n           [--message-thresholds xml_file]"
   << "\n           [--unphysical-event-mask mask]"
   << "\n           [--event-record-print-level level]"
   << "\n           [--mc-job-status-refresh-rate  rate]"
   << "\n           [--cache-file root_file]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n";
}
//________________________________________________________________________________________
