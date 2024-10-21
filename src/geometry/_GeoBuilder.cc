// G4 material
#include "G4Material.hh"
#include "G4NistManager.hh"
// G4 geometry
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
// G4 visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
// G4 constants/units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// User include:
#include "geometry/_GeoBuilder.hh"

namespace MuGeoBuilder
{
    // ----------------------------------------------------------------------
    // Set Construction Materials, saved in namespace Material
    const auto _nist_ = G4NistManager::Instance();
    G4Element *Material::H = _nist_->FindOrBuildElement("H");
    G4Element *Material::C = _nist_->FindOrBuildElement("C");
    G4Element *Material::N = _nist_->FindOrBuildElement("N");
    G4Element *Material::O = _nist_->FindOrBuildElement("O");
    G4Element *Material::F = _nist_->FindOrBuildElement("F");
    G4Element *Material::Al = _nist_->FindOrBuildElement("Al");
    G4Element *Material::Si = _nist_->FindOrBuildElement("Si");
    G4Element *Material::S = _nist_->FindOrBuildElement("S");
    G4Element *Material::Ar = _nist_->FindOrBuildElement("Ar");
    G4Element *Material::Ca = _nist_->FindOrBuildElement("Ca");
    G4Material *Material::Air = _nist_->FindOrBuildMaterial("G4_AIR");
    G4Material *Material::Aluminum = _nist_->FindOrBuildMaterial("G4_Al");
    G4Material *Material::Bakelite = _nist_->FindOrBuildMaterial("G4_BAKELITE");
    G4Material *Material::Copper = _nist_->FindOrBuildMaterial("G4_Cu");
    G4Material *Material::Concrete = _nist_->FindOrBuildMaterial("G4_CONCRETE");
    G4Material *Material::Iron = _nist_->FindOrBuildMaterial("G4_Fe");
    G4Material *Material::PolystyreneFoam = _nist_->BuildMaterialWithNewDensity("PolystyreneFoam", "G4_POLYSTYRENE", 32.0 * kg / m3);
    G4Material *Material::Polyvinyltoluene = _nist_->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material *Material::LiquidArgon = new G4Material("LiquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
    G4Material *Material::Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);// Vacuum


    // ----------------------------------------------------------------------
    // Geometry Builder Class
    // Set messenger directory
    const std::string Builder::MessengerDirectory = "/det/";

    // Constructor
    Builder::Builder(const std::string &detector_name) : G4UImessenger(MessengerDirectory + detector_name, "Detector detector_name")
    {

        // Make a few more material. This couldn't be down outside the class
        Material::PlasticScintillator = new G4Material("PlasticScintillator", 1.032*g/cm3, 2);
        Material::PlasticScintillator->AddElement(Material::C, 9);
        Material::PlasticScintillator->AddElement(Material::H, 10);    
    }

    G4VPhysicalVolume *Builder::Construct()
    {
    }

    // Messenger related
    void Builder::SetNewValue(G4UIcommand *command, G4String value)
    {
        (void)(command);
        (void)(value);
    }


} // namespace MuGeoBuilder