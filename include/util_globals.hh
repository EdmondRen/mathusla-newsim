#ifndef UTIL_GLOBALS_HH
#define UTIL_GLOBALS_HH

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
#include <variant>
#include <map>
#include <iomanip>

namespace util
{
    static std::string VERSION = "1.0";

    namespace globals
    {
        extern std::string PROJECT_SOURCE_DIR;

    } // namespace globals
} // namespace util

using util::globals::PROJECT_SOURCE_DIR;

#endif