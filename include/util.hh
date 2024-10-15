// util.h
#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
#include <filesystem>

namespace util
{

    namespace path
    {
        std::filesystem::path getExecutablePath();
    } // namespace path    
}

#endif // UTIL_H