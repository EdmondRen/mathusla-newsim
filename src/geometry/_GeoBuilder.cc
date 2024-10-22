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
    // Geometry Builder Class
    Builder::Builder() {}
    G4VPhysicalVolume *Builder::Construct() { return 0; }

    // ----------------------------------------------------------------------
    // Helper functions for building geometry, such as
    //   * materials
    //   * solids
    //   * visualization styles
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
    G4Material *Material::LiquidArgon = new G4Material("LiquidArgon", 18., 39.95 * g / mole, 1.390 * g / cm3);
    G4Material *Material::Vacuum = new G4Material("Galactic", 1., 1.01 * g / mole, universe_mean_density, kStateGas, 2.73 * kelvin, 3.e-18 * pascal); // Vacuum
    G4Material *Material::PlasticScintillator = new G4Material("PlasticScintillator", 1.032 * g / cm3, 2);
    G4Material *Material::CaCO3 = new G4Material("CaCO3", 2.71*g/cm3, 3);
    G4Material *Material::Kaolinite = new G4Material("Clay", 2.65*g/cm3, 4);
    G4Material *Material::SiO2 = new G4Material("Quartz", 2.445*g/cm3, 2);
    G4Material *Material::Marl = new G4Material("Marl", 2.46*g/cm3, 2);
    G4Material *Material::GroundMix = new G4Material("GroundMix", 2.54*g/cm3, 2);
    // Set the details of thoses composite material. This couldn't be down outside of a function
    G4Material *make_materials()
    {
        Material::PlasticScintillator->AddElement(Material::C, 9);
        Material::PlasticScintillator->AddElement(Material::H, 10);
        Material::CaCO3->AddElement(Material::Ca, 1);
        Material::CaCO3->AddElement(Material::C,  1);
        Material::CaCO3->AddElement(Material::O,  3);
        Material::Kaolinite->AddElement(Material::Al, 2);
        Material::Kaolinite->AddElement(Material::Si, 2);
        Material::Kaolinite->AddElement(Material::O,  9);
        Material::Kaolinite->AddElement(Material::H,  4);
        Material::SiO2->AddElement(Material::Si, 1);
        Material::SiO2->AddElement(Material::O, 2);
        Material::Marl->AddMaterial(Material::Kaolinite, 35*perCent);
        Material::Marl->AddMaterial(Material::CaCO3,     65*perCent);
        Material::GroundMix->AddMaterial(Material::Marl, 50*perCent);
        Material::GroundMix->AddMaterial(Material::SiO2, 50*perCent);  
        return Material::PlasticScintillator;     
    }
    G4Material *temp = make_materials();

    namespace Vis
    {
        //__Sensitive Material Attribute Definition_____________________________________________________
        const G4VisAttributes SensitiveAttributes1()
        {
            auto attr = G4VisAttributes(G4Colour(0., 1., 0., 1.0));
            attr.SetForceSolid(true);
            return attr;
        }

        const G4VisAttributes SensitiveAttributes2()
        {
            auto attr = G4VisAttributes(G4Colour(1., 0., 0., 1.0));
            attr.SetForceSolid(true);
            return attr;
        }

        const G4VisAttributes HighlightRed()
        {
            auto attr = G4VisAttributes(G4Colour(1., 0., 0.));
            attr.SetForceSolid(true);
            return attr;
        }

        const G4VisAttributes TransparentBrown()
        {
            auto attr = G4VisAttributes(G4Colour(0.8, 0.5,0.4, 0.3));
            attr.SetForceSolid(true);
            return attr;
        }        


        //__Casing Material Attribute Definition________________________________________________________
        const G4VisAttributes CasingAttributes()
        {
            auto attr = G4VisAttributes(G4Colour(0., 0., 1., 0.2));
            attr.SetForceWireframe(true);
            return attr;
        }

        //__Border Attribute Definition_________________________________________________________________
        const G4VisAttributes BorderAttributes()
        {
            return G4VisAttributes(false);
        }

        //__Iron Plating Attributes
        const G4VisAttributes IronAttributes()
        {
            auto attr = G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
            attr.SetForceSolid(true);
            return attr;
        }

        const G4VisAttributes AlAttributes()
        {
            auto attr = G4VisAttributes(G4Colour(0.4, 0.6, 0.7));
            attr.SetForceSolid(true);
            return attr;
        }


    } // namespace Vis

} // namespace MuGeoBuilder