#include "generator/_Generator.hh"

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

    // Other helper functions
    std::ostream &Generator::Print(std::ostream &os) const { return os; }

}