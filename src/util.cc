#include "util.hh"

#ifdef _BSD_SOURCE
#include <sysexits.h>
#else
#define EX_USAGE EXIT_FAILURE
#endif

#include <vector>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <limits.h>
#include <unordered_map>
#include <memory>
#include <stdexcept>

namespace util
{

    namespace py
    {

    } // namespace py

    namespace path
    {
        /* Get the path of the current executable */
        std::filesystem::path getExecutablePath()
        {
            char buffer[PATH_MAX];
            ssize_t count = readlink("/proc/self/exe", buffer, PATH_MAX);
            if (count != -1)
            {
                buffer[count] = '\0';
                return std::filesystem::path(buffer);
            }
            return std::filesystem::path();
        }
    } // namespace path

    namespace io
    {
        std::string readFileToString(const std::string &filename)
        {
            std::ifstream file(filename);
            if (!file.is_open())
            {
                throw std::runtime_error("Could not open file: " + filename);
            }

            std::ostringstream oss;
            oss << file.rdbuf(); // Read the file's buffer into the string stream

            return oss.str(); // Return the contents as a string
        }

        std::string readFileToString_CRY(const std::string &filename)
        {
            // Read the cry input file
            std::ifstream inputFile(filename);
            char buffer[1000];

            std::string setupString("");

            if (inputFile.fail())
            {
                throw std::runtime_error("Could not open file: " + filename);
            }
            else
            {
                while (!inputFile.getline(buffer, 1000).eof())
                {
                    setupString.append(buffer);
                    setupString.append(" ");
                }
            }

            return setupString;
        }

    } // namespace io

} // namespace util
