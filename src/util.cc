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

} // namespace util
