#ifndef GeoBuilder_HH
#define GeoBuilder_HH 1

#include "G4Material.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GenericMessenger.hh"
#include <G4UIcommand.hh>
#include <G4VisAttributes.hh>

namespace MuGeoBuilder
{

    // "y_side_direction": int (one of 0,1,2)
    // "z_side_direction": int (one of 0,1,2)
    // "y_side_coord": float
    // "z_side_coord": float
    struct BarPosition
    {
        G4ThreeVector y_side_direction;
        G4ThreeVector z_side_direction;
        G4ThreeVector bar_center_coord;

        BarPosition(G4ThreeVector _y_side_direction,
                    G4ThreeVector _z_side_direction,
                    G4ThreeVector _bar_center_coord): y_side_direction(_y_side_direction.x(), _y_side_direction.y(), _y_side_direction.z()),
                    z_side_direction(_z_side_direction.x(), _z_side_direction.y(), _z_side_direction.z()),
                    bar_center_coord(_bar_center_coord.x(), _bar_center_coord.y(), _bar_center_coord.z())
        {
        }

        // BarPosition(G4ThreeVector _y_side_direction,
        //             G4ThreeVector _z_side_direction,
        //             G4ThreeVector _bar_center_coord)
        // {
        //     y_side_direction = _y_side_direction;
        //     z_side_direction = _z_side_direction;
        //     bar_center_coord = _bar_center_coord;
        // }        
    };

    // Geometry Builder Class
    // This class takes care of detector construction
    // as well as assigning sensitive detector to correct volumes
    class Builder
    {
    public:
        // Builder();
        // virtual ~Builder() =0;

        // Core function 1:
        // Construct physics volume, should return world PV
        virtual G4VPhysicalVolume *Construct() = 0;
        G4VPhysicalVolume *worldPV;

        // Core function 2:
        // Set the sensitive detector for this geometry
        virtual void ConstructSD(G4VSensitiveDetector *detector) = 0;
        G4VSensitiveDetector *fdetector;

        // Core function 3:
        // (For digitizer) Get a unique detector ID for each bar based on the copy number
        // Depends on the implementation of the geometry, local coordinates may be needed.
        // For example, if the bottom-level touchable is a monolithic plate instead of individual bars...
        virtual long long int GetDetectorID(std::vector<int> copy_numbers, G4ThreeVector local_coords) = 0;

        // Core function 4:
        // (For digitizer) Make a map from detector ID to bar information dict
        virtual BarPosition GetBarPosition(long long detector_id) = 0;

        // Utilities
        // Find the position of a physical volume in world
        std::vector<CLHEP::Hep3Vector> GetPositionInWorld(const G4String &physName);

        // Varibles
        std::string DetectorName;
        // Messenger related
        G4GenericMessenger *fMessenger;

    private:
    };

    namespace Material
    {
        extern G4Element *H;
        extern G4Element *C;
        extern G4Element *N;
        extern G4Element *O;
        extern G4Element *F;
        extern G4Element *Al;
        extern G4Element *Si;
        extern G4Element *S;
        extern G4Element *Ar;
        extern G4Element *Ca;
        extern G4Material *Air;
        extern G4Material *Aluminum;
        extern G4Material *Bakelite;
        extern G4Material *Copper;
        extern G4Material *Concrete;
        extern G4Material *Iron;
        extern G4Material *PolystyreneFoam;
        extern G4Material *Polyvinyltoluene;
        extern G4Material *PlasticScintillator;
        extern G4Material *LiquidArgon;
        extern G4Material *Vacuum;
        extern G4Material *CaCO3;
        extern G4Material *Kaolinite;
        extern G4Material *SiO2;
        extern G4Material *Marl;
        extern G4Material *GroundMix;

    } /* namespace Material */

    namespace Vis
    {
        extern std::map<std::string, G4VisAttributes> styles;

    } /* namespace Vis*/

} // namespace MuGeoBuilder

#endif
