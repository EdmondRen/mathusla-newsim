#include "generator/_Generator.hh"
#include "G4RandomTools.hh" // For G4UniformRand

namespace MuGenerators
{

    //__Generator Messenger Directory Path__________________________________________________________
    const std::string Generator::MessengerDirectory = "/gen/";

    Generator::Generator(const std::string &name,
                         const std::string &description) : G4UImessenger(MessengerDirectory + name, description)
    {
    }

    void Generator::GeneratePrimaryVertex(G4Event *anEvent) { (void)(anEvent); }

    // Core function 2: GeneratePrimaryVertex()
    // This is used to set generator parameters
    void Generator::SetNewValue(G4UIcommand *command,
                                G4String value)
    {
        (void)(command);
        (void)(value);
    }

    //--------------------------------------------------------------------------
    // Other helper functions
    std::ostream &Generator::Print(std::ostream &os) const { return os; }

    int Generator::GetEntries() const
    {
        return -1;
    }

    // Wrapper for random number generator
    double GenerateRandomInRange(double min, double max)
    {
        if (min >= max)
        {
            throw std::invalid_argument("Invalid range: min must be less than max");
        }
        // Generate a random number in the range [min, max)
        return min + (max - min) * G4UniformRand();
    }

    // Helper function to find intersection range
    bool intersectSlab(double p0, double d, double min, double max, double &tmin, double &tmax)
    {
        if (std::abs(d) < 1e-8)
        {
            // Line is parallel to the slab
            return p0 >= min && p0 <= max;
        }
        double t1 = (min - p0) / d;
        double t2 = (max - p0) / d;
        if (t1 > t2)
            std::swap(t1, t2);
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        return tmin <= tmax;
    }

    bool doesLineIntersectBox(const Vec3 &p0, const Vec3 &d, const Vec3 &boxMin, const Vec3 &boxMax)
    {
        double tmin = -INFINITY, tmax = INFINITY;

        // Check x-axis slab
        if (!intersectSlab(p0.x, d.x, boxMin.x, boxMax.x, tmin, tmax))
            return false;

        // Check y-axis slab
        if (!intersectSlab(p0.y, d.y, boxMin.y, boxMax.y, tmin, tmax))
            return false;

        // Check z-axis slab
        if (!intersectSlab(p0.z, d.z, boxMin.z, boxMax.z, tmin, tmax))
            return false;

        return true;
    }    
    //--------------------------------------------------------------------------

}