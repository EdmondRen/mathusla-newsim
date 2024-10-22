// util.h
#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
#include <filesystem>
#include <unordered_map>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <any>

namespace util
{   
    namespace globals
    {
        extern std::string PROJECT_SOURCE_DIR;

    } // namespace globals

    /*
    Python-like interface, such as dict, print...
    */
    namespace py
    {
        // Python-like dictionary class
        /*
        Example:
        #include "util.hh"
        int main() {
            util::py::dict my_dict;
            // Inserting values of different types
            my_dict.set("int_key", 42);
            my_dict.set("double_key", 3.14159);
            my_dict.set("string_key", std::string("Hello, World"));
            my_dict.set("vector_key", std::vector<int>{1, 2, 3, 4, 5});
            // Printing the values
            auto int_key = my_dict.get("int_key", int())
            std::cout << int_key << std::endl;
            return 0;
        }
        */
        class dict {
        public:
            std::unordered_map<std::string, std::any> _tVars;

            template <typename T>
            T get(const std::string &key, T defaultValue) const
            {
                auto it = _tVars.find(key);
                if (it == _tVars.end())
                    return defaultValue;

                return std::any_cast<T>(it->second);
            }

            template <typename T>
            void set(const std::string &key, T value)
            {
                _tVars[key] = value;
            }
        };       

    } // namespace py

    /*
    Path-related helper function
    */
    namespace path
    {
        std::filesystem::path getExecutablePath();
    } // namespace path

    namespace io
    {
        std::string readFileToString(const std::string &filename);
    } // namespace io
    
}

#endif // UTIL_H




//   util::py::argparse parser;
//   // Add positional and optional arguments
//   parser.addOptional<int>("--threads", 1);
//   parser.addOptional<G4String>("--session", "MathuslaSim");
//   parser.addOptional<std::string>("--macro", "", true); // True: Multiple values for "--macro"

//   util::py::dict args;
//   parser.parse(argc, argv, args);

//   // * macro: need to separate the macro name and forwarding arguments
//   std::vector<std::string> macro_commands;
//   std::string macro;
//   macro_commands = args.contains("--macro") ? args.get<std::vector<std::string> >("--macro") : macro_commands;
//   if (macro_commands.size()){
//     macro = macro_commands[0];
//   }
//   auto session = args.contains("--session") ? args.get<G4String>("session") : "";
//   auto nThreads = args.contains("--threads") ? args.get<G4int>("threads") : 1;