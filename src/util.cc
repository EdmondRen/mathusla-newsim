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

#include <cstdio>
#include <sys/stat.h>
#if defined(_WIN32)
#include <windows.h>
#endif

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

        //__Create Directory___
        int create_directory(const std::string &newDirPath)
        {
#if defined(_WIN32)
            return _mkdir(dir.c_str());
#else
            // Get current directory. If it is relative, use the  directory of the executable
            auto testDir = newDirPath[0] == '/' ? newDirPath : path::getExecutablePath().string();
            // So that we can make new directory with same permissions
            auto path_current = std::filesystem::path(testDir);
            auto path_parent = path_current.parent_path().string();
            // Get the parent's status
            struct stat parentStat;
            if (stat(path_parent.c_str(), &parentStat) != 0)
            {
                std::cout << "Error getting parent directory status: " << strerror(errno) << std::endl;
                return false;
            }
            // Extract the parent's permissions
            mode_t parentPermissions = parentStat.st_mode & (S_IRWXU | S_IRWXG | S_IRWXO); // Read, write, execute for user, group, others

            // Finally, create the new directory with the same permissions as the parent
            if (mkdir(newDirPath.c_str(), parentPermissions) != 0)
            {
                std::cout << "Error creating directory: " << strerror(errno) << std::endl;
                return false;
            }
            return true;

#endif
        }

        ParHandler::ParHandler(std::string par_card_str) : filename(par_card_str)
        {

            std::string par_path = std::string(__FILE__);

            std::ifstream infile(par_card_str);

            if (infile.is_open())
            {
                file_opened = true;

                std::string line;
                while (std::getline(infile, line))
                {
                    if (line[0] == '#')
                        continue;

                    std::istringstream iss(line);
                    std::string par;
                    double value;

                    if (!(iss >> par >> value))
                    {
                        continue;
                    } // error

                    config[par] = value;
                }
            }
            else
            {
                file_opened = false;
                std::cout << "Error opening parameter card file" << std::endl;
            }
        }

        std::string ParHandler::GetString()
        {
            // Create an input file stream
            std::ifstream file(filename);

            // Check if the file was opened successfully
            if (!file.is_open())
            {
                throw std::runtime_error("Could not open file: " + filename);
            }

            // Use a stringstream to read the file contents
            std::ostringstream sstr;
            sstr << file.rdbuf();

            // Return the contents as a string
            return sstr.str();
        }

    } // namespace io

    namespace notstd
    {
        // a method to check if a substring is in a string.
        const in_helper::in_t in{}; // Define `in` here

        bool strfind(const std::string &substring, const std::string &str)
        {
            // Use std::string::find to check if the substring exists
            std::size_t found = str.find(substring);
            return found != std::string::npos;
        }
    }

    namespace vector
    {
        std::vector<std::vector<int>> splitVectorByDelimiter(const std::vector<int> &input, int delimiter)
        {
            std::vector<std::vector<int>> result;
            std::vector<int> currentSubVector;

            for (int value : input)
            {
                if (value == delimiter)
                {
                    // If we encounter the delimiter, push the current sub-vector to the result.
                    if (!currentSubVector.empty())
                    {
                        result.push_back(currentSubVector);
                        currentSubVector.clear();
                    }
                }
                else
                {
                    // Otherwise, add the value to the current sub-vector.
                    currentSubVector.push_back(value);
                }
            }

            // If there's anything left in the current sub-vector, push it to the result.
            if (!currentSubVector.empty())
            {
                result.push_back(currentSubVector);
            }

            return result;
        }
    }

} // namespace util

// using  util::notstd::in;
// using  util::py::print;
