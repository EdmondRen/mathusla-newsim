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

#include "G4TransportationManager.hh"
#include "G4PhysicalVolumesSearchScene.hh"

// User include:
#include "geometry/_GeoBuilder.hh"

namespace MuGeoBuilder
{

    // ----------------------------------------------------------------------
    // Geometry Builder Class
    // Builder::Builder() {}
    // G4VPhysicalVolume *Builder::Construct() { return 0; }
    // void Builder::ConstructSD(G4VSensitiveDetector *detector) { (void)detector; }

    std::vector<CLHEP::Hep3Vector> Builder::GetPositionInWorld(const G4String &physName)
    {
        G4TransportationManager *transportationManager =
            G4TransportationManager::GetTransportationManager();

        size_t nWorlds = transportationManager->GetNoWorlds();

        std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;
        std::vector<G4VPhysicalVolume *>::iterator iterWorld =
            transportationManager->GetWorldsIterator();
        for (size_t i = 0; i < nWorlds; ++i, ++iterWorld)
        {
            G4PhysicalVolumeModel searchModel(*iterWorld); // Unlimited depth.
            G4ModelingParameters mp;                       // Default - no culling.
            searchModel.SetModelingParameters(&mp);
            G4PhysicalVolumesSearchScene searchScene(&searchModel, physName);
            searchModel.DescribeYourselfTo(searchScene); // Initiate search.

            for (const auto &findings : searchScene.GetFindings())
            {
                findingsVector.push_back(findings);
            }
        }

        for (const auto &findings : findingsVector)
        {
            G4cout
                << findings.fFoundBasePVPath
                << " " << findings.fpFoundPV->GetName()
                << " " << findings.fFoundPVCopyNo
                << " (mother logical volume: "
                << findings.fpFoundPV->GetMotherLogical()->GetName()
                << ")"
                << G4endl;
        }

        std::vector<CLHEP::Hep3Vector> vPos;

        if (!findingsVector.size())
        {
            G4cout << physName << " not found" << G4endl;
        }
        else
        {
            for (const auto &findings : findingsVector)
            {
                vPos.push_back(findings.fFoundObjectTransformation.getTranslation());
            }

            for (const auto &pos : vPos)
            {
                G4cout << pos / m << G4endl;
            }
        }

        return vPos;
    }

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
    G4Material *Material::CaCO3 = new G4Material("CaCO3", 2.71 * g / cm3, 3);
    G4Material *Material::Kaolinite = new G4Material("Clay", 2.65 * g / cm3, 4);
    G4Material *Material::SiO2 = new G4Material("Quartz", 2.445 * g / cm3, 2);
    G4Material *Material::Marl = new G4Material("Marl", 2.46 * g / cm3, 2);
    G4Material *Material::GroundMix = new G4Material("GroundMix", 2.54 * g / cm3, 2);
    // Set the details of thoses composite material. This couldn't be down outside of a function
    G4Material *make_materials()
    {
        Material::PlasticScintillator->AddElement(Material::C, 9);
        Material::PlasticScintillator->AddElement(Material::H, 10);
        Material::CaCO3->AddElement(Material::Ca, 1);
        Material::CaCO3->AddElement(Material::C, 1);
        Material::CaCO3->AddElement(Material::O, 3);
        Material::Kaolinite->AddElement(Material::Al, 2);
        Material::Kaolinite->AddElement(Material::Si, 2);
        Material::Kaolinite->AddElement(Material::O, 9);
        Material::Kaolinite->AddElement(Material::H, 4);
        Material::SiO2->AddElement(Material::Si, 1);
        Material::SiO2->AddElement(Material::O, 2);
        Material::Marl->AddMaterial(Material::Kaolinite, 35 * perCent);
        Material::Marl->AddMaterial(Material::CaCO3, 65 * perCent);
        Material::GroundMix->AddMaterial(Material::Marl, 50 * perCent);
        Material::GroundMix->AddMaterial(Material::SiO2, 50 * perCent);
        return Material::PlasticScintillator;
    }
    G4Material *temp = make_materials();

    namespace Vis
    {
        // make visualization styles
        std::map<std::string, G4VisAttributes> make_styles()
        {
            std::map<std::string, G4VisAttributes> styles_store;

            styles_store["SensitiveAttributes1"] = G4VisAttributes(G4Colour(0.3, 1., 0., 0.9));
            styles_store["SensitiveAttributes1"].SetForceSolid(true);

            styles_store["SensitiveAttributes2"] = G4VisAttributes(G4Colour(0.1, 0.5, 0., 1));
            styles_store["SensitiveAttributes2"].SetForceSolid(true);

            styles_store["HighlightRed"] = G4VisAttributes(G4Colour(1, 0, 0, 1));
            styles_store["HighlightRed"].SetForceSolid(true);

            styles_store["TransparentBrown"] = G4VisAttributes(G4Colour(0.8, 0.5, 0.4, 0.3));
            styles_store["TransparentBrown"].SetForceSolid(true);

            styles_store["TransparentBlue"] = G4VisAttributes(G4Colour(0.2, 0.6, 1, 0.1));
            styles_store["TransparentBlue"].SetForceSolid(true);

            styles_store["CasingAttributes"] = G4VisAttributes(G4Colour(0., 0., 1., 0.2));
            styles_store["CasingAttributes"].SetForceWireframe(true);

            styles_store["SensitiveAttributes_border"] = G4VisAttributes(G4Colour(0.1, 0.5, 0.0));
            styles_store["SensitiveAttributes_border"].SetVisibility(true);
            styles_store["SensitiveAttributes_border"].SetForceWireframe(true);
            styles_store["SensitiveAttributes_border"].SetLineWidth(3);

            styles_store["IronAttributes"] = G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 1.0));
            styles_store["IronAttributes"].SetForceSolid(true);

            styles_store["AlAttributes"] = G4VisAttributes(G4Colour(0.4, 0.6, 0.7, 0.4));
            styles_store["AlAttributes"].SetForceSolid(true);

            // A few general styles
            styles_store["Invisible"] = G4VisAttributes(G4Colour(1.0, 1.0, 1.0)); // White color for clarity
            styles_store["Invisible"].SetVisibility(false);

            styles_store["WhiteFrame"] = G4VisAttributes(G4Colour(0., 0., 1., 1.));
            styles_store["WhiteFrame"].SetForceWireframe(true);     

            styles_store["WhiteCloud"] = G4VisAttributes(G4Colour(0., 0., 1.));
            styles_store["WhiteCloud"].SetVisibility(true);     
            styles_store["WhiteCloud"].SetForceSolid(false);     
            styles_store["WhiteCloud"].SetForceWireframe(true);     
            // styles_store["WhiteCloud"].SetForcePoints(true);     

            return styles_store;
        }

        std::map<std::string, G4VisAttributes> styles = make_styles();

    } // namespace Vis

} // namespace MuGeoBuilder