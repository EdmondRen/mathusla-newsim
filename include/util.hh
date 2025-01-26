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
#include <variant>
#include <map>
#include <iomanip>

namespace util
{
    static std::string VERSION = "0.1";

    namespace globals
    {
        extern std::string PROJECT_SOURCE_DIR;

    } // namespace globals

    /*
    Python-like interface, such as dict, print...
    */
    namespace py
    {

        // Base template for variadic arguments
        template <typename T>
        void print_item(std::ostringstream &oss, const T &item)
        {
            oss << item;
        }

        // Overload for vector printing
        template <typename T>
        void print_item(std::ostringstream &oss, const std::vector<T> &vec)
        {
            oss << "[";
            for (size_t i = 0; i < vec.size(); ++i)
            {
                oss << vec[i];
                if (i != vec.size() - 1)
                {
                    oss << ", ";
                }
            }
            oss << "]";
        }

        // python-like print function
        //  - Usage: print(arg1, arg2, ...)
        //  - Does not support different sep and end
        //    Variadic template function to handle any number of arguments
        template <typename... Args>
        void print(Args... args)
        {
            // Default values for sep and end
            std::string sep = " ";
            std::string end = "\n";

            std::ostringstream oss;

            // Using fold expression to unpack the arguments and insert the separator
            ((print_item(oss, args), oss << sep), ...);

            std::string output = oss.str();

            // Remove the last separator
            if (!output.empty())
            {
                output.erase(output.size() - sep.size());
            }

            // Print the result with the specified ending
            std::cout << output << end;
        }

        template <typename... Args>
        std::string fprint(Args... args)
        {
            // Default values for sep and end
            std::string sep = " ";
            std::string end = "\n";

            std::ostringstream oss;

            // Using fold expression to unpack the arguments and insert the separator
            ((print_item(oss, args), oss << sep), ...);

            std::string output = oss.str();

            // Remove the last separator
            if (!output.empty())
            {
                output.erase(output.size() - sep.size());
            }

            // Print the result with the specified ending
            output += end;
            return output;
        }

        // Base case for variadic template: no arguments left to process
        inline void format_impl(std::ostringstream &oss, const std::string &format)
        {
            oss << format;
        }

        // Helper function to apply formatting
        template <typename T>
        void apply_format(std::ostringstream &oss, const std::string &spec, const T &value)
        {
            if (spec.empty())
            {
                oss << value; // No formatting specified
                return;
            }

            int width = 0;      // Default width
            int precision = -1; // Default precision (unset)
            size_t dot_pos = std::string::npos;

            // Skip leading ':' if present
            std::string clean_spec = spec[0] == ':' ? spec.substr(1) : spec;

            try
            {
                dot_pos = clean_spec.find('.');

                if (dot_pos != std::string::npos)
                {
                    // Parse width and precision separately
                    if (dot_pos > 0)
                    {
                        width = std::stoi(clean_spec.substr(0, dot_pos));
                    }
                    precision = std::stoi(clean_spec.substr(dot_pos + 1));
                }
                else if (!clean_spec.empty())
                {
                    // Only width specified
                    width = std::stoi(clean_spec);
                }
            }
            catch (const std::invalid_argument &)
            {
                throw std::invalid_argument("Invalid format specifier: " + spec);
            }
            catch (const std::out_of_range &)
            {
                throw std::out_of_range("Format specifier out of range: " + spec);
            }

            // Apply width and precision
            if (precision >= 0)
            {
                oss << std::fixed << std::setprecision(precision);
            }
            if (width > 0)
            {
                oss << std::setw(width) << std::setfill(' ');
            }

            oss << value;
        }

        // Variadic template to process each argument
        template <typename T, typename... Args>
        void format_impl(std::ostringstream &oss, const std::string &format, const T &value, const Args &...args)
        {
            size_t open_brace = format.find('{');
            size_t close_brace = format.find('}', open_brace);

            if (open_brace == std::string::npos || close_brace == std::string::npos || close_brace < open_brace)
            {
                throw std::invalid_argument("Invalid format string: unmatched braces or missing arguments");
            }

            // Append the text before the placeholder
            oss << format.substr(0, open_brace);

            // Extract the specifier inside the braces (e.g., `.2f`)
            std::string spec = format.substr(open_brace + 1, close_brace - open_brace - 1);

            // Format and insert the value
            apply_format(oss, spec, value);

            // Process the remaining format string after the placeholder
            format_impl(oss, format.substr(close_brace + 1), args...);
        }

        // Python-like format string
        // Main format function
        template <typename... Args>
        std::string f(const std::string &format_str, const Args &...args)
        {
            std::ostringstream oss;
            format_impl(oss, format_str, args...);
            return oss.str();
        }

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
        enum KEYTYPE
        {
            __single__,
            __vector__,
        };

        // Tuple can hold data types of I (int), F(float), D(double) or S(string)
        enum DTYPE
        {
            __int__,
            __float__,
            __double__,
            __string__,
        };

        // Define a variant type that can hold either an int or a std::string
        using DictVariant = std::variant<int,
                                         float,
                                         double,
                                         std::string,
                                         std::vector<int> *,
                                         std::vector<float> *,
                                         std::vector<double> *>;

        class Column
        {
        public:
            // Constructor
            Column(std::string kname, int kindex, KEYTYPE ktype, DTYPE dtype)
            {
                this->key_type = ktype;
                this->data_type = dtype;
                this->key_index = kindex;
                this->key_name = kname;
            }

            // Destructor
            ~Column()
            {
            }
            std::string key_name;
            int key_index;
            enum KEYTYPE key_type;
            enum DTYPE data_type;

            int data_int;
            float data_float;
            double data_double;
            std::string data_string;
            std::vector<int> data_vec_int;
            std::vector<float> data_vec_float;
            std::vector<double> data_vec_double;

            void push_back(const double value)
            {
                if (this->data_type == __int__)
                    data_vec_int.push_back((int)value);
                else if (this->data_type == __float__)
                    data_vec_float.push_back((float)value);
                else if (this->data_type == __double__)
                {
                    data_vec_double.push_back(value);
                }
            }

            void clear()
            {
                if (this->data_type == __int__)
                    data_vec_int.clear();
                else if (this->data_type == __float__)
                    data_vec_float.clear();
                else if (this->data_type == __double__)
                    data_vec_double.clear();

                data_int = 0;
                data_float = 0;
                data_double = 0;
                data_string = "";
            }

            // Overload assignment operator for MyClass
            Column &operator=(const Column &other)
            {
                if (this != &other)
                { // Self-assignment check
                    data_int = other.data_int;
                    data_float = other.data_float;
                    data_double = other.data_double;
                }
                return *this;
            }

            Column &operator=(const double value)
            {
                if (this->data_type == __int__)
                    data_int = (int)value;
                else if (this->data_type == __float__)
                    data_float = (float)value;
                else if (this->data_type == __double__)
                    data_double = value;
                return *this;
            }

            Column &operator=(const std::string value)
            {
                data_string = value;
                return *this;
            }
        };

        class Dict
        {
        public:
            using data_t = std::map<std::string, Column *>;
            data_t data;

            ~Dict()
            {
            }

            void clear()
            {
                data_t::iterator it = this->data.begin();
                while (it != this->data.end())
                {
                    it->second->clear();
                    it++;
                }
            }

            std::vector<std::string> keys()
            {
                std::vector<std::string> keys;
                data_t::iterator it = this->data.begin();
                while (it != this->data.end())
                {
                    keys.push_back(it->first);
                    ++it;
                }
                return keys;
            }

            int size()
            {
                return ncolumns;
            }

            // Overload the subscript operator for const objects
            Column &operator[](const std::string &key_name)
            {
                return Get(key_name);
            }

            Column &operator[](int key_ind)
            {
                return Get(key_ind);
            }

            Column &Get(std::string key_name)
            {
                return *this->data[key_name];
            }

            Column &Get(int key_ind)
            {
                return *this->data[keys_map[key_ind]];
            }

            void Add(std::string key_name, KEYTYPE key_type, DTYPE data_type)
            {
                this->data[key_name] = new Column(key_name, ncolumns, key_type, data_type);
                this->keys_map.push_back(key_name);
                ncolumns += 1;
            }

            // clang-format off
            // DictVariant Get(std::string key)
            // {
            //     if (this->data[key]->key_type == __single__ && this->data[key]->data_type == __int__)
            //         return this->data[key]->data_int;
            //     else if (this->data[key]->key_type == __single__ && this->data[key]->data_type == __float__)
            //         return this->data[key]->data_float;
            //     else if (this->data[key]->key_type == __single__ && this->data[key]->data_type == __double__)
            //         return this->data[key]->data_double;
            //     else if (this->data[key]->key_type == __single__ && this->data[key]->data_type == __string__)
            //         return this->data[key]->data_double;
            //     else if (this->data[key]->key_type == __vector__ && this->data[key]->data_type == __int__)
            //         return & this->data[key]->data_vec_int;
            //     else if (this->data[key]->key_type == __vector__ && this->data[key]->data_type == __float__)
            //         return & this->data[key]->data_vec_float;
            //     else if (this->data[key]->key_type == __vector__ && this->data[key]->data_type == __double__)
            //         return & this->data[key]->data_vec_double;
            //     return 999999;
            // }
            // clang-format on

        private:
            int ncolumns;
            std::vector<std::string> keys_map;
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
        std::string readFileToString_CRY(const std::string &filename);

        // Create Directory with same permission as its parent_
        int create_directory(const std::string &newDirPath);

        class ParHandler
        {

        public:
            std::map<std::string, double> config;
            bool file_opened;
            ParHandler(std::string filename);

            std::map<std::string, double> &GetConfig() { return config; };

            double &operator[](const std::string &key_name)
            {
                return this->config[key_name];
            }
        };

    } // namespace io

    namespace notstd
    {
        namespace ca_helper
        {
            template <template <class...> class, class, class...>
            struct can_apply : std::false_type
            {
            };
            template <class...>
            struct voider
            {
                using type = void;
            };
            template <class... Ts>
            using void_t = typename voider<Ts...>::type;

            template <template <class...> class Z, class... Ts>
            struct can_apply<Z, void_t<Z<Ts...>>, Ts...> : std::true_type
            {
            };
        }
        template <template <class...> class Z, class... Ts>
        using can_apply = ca_helper::can_apply<Z, void, Ts...>;

        namespace find_helper
        {
            template <class C, class T>
            using dot_find_r = decltype(std::declval<C>().find(std::declval<T>()));
            template <class C, class T>
            using can_dot_find = can_apply<dot_find_r, C, T>;
            template <class C, class T>
            constexpr std::enable_if_t<can_dot_find<C &, T>{}, bool>
            find(C &&c, T &&t)
            {
                using std::end;
                return c.find(std::forward<T>(t)) != end(c);
            }
            template <class C, class T>
            constexpr std::enable_if_t<!can_dot_find<C &, T>{}, bool>
            find(C &&c, T &&t)
            {
                using std::begin;
                using std::end;
                return std::find(begin(c), end(c), std::forward<T>(t)) != end(c);
            }
            template <class C, class T>
            constexpr bool finder(C &&c, T &&t)
            {
                return find(std::forward<C>(c), std::forward<T>(t));
            }
        }
        template <class C, class T>
        constexpr bool find(C &&c, T &&t)
        {
            return find_helper::finder(std::forward<C>(c), std::forward<T>(t));
        }
        struct finder_t
        {
            template <class C, class T>
            constexpr bool operator()(C &&c, T &&t) const
            {
                return find(std::forward<C>(c), std::forward<T>(t));
            }
            constexpr finder_t() {}
        };
        constexpr finder_t finder{};
        namespace named_operator
        {
            template <class D>
            struct make_operator
            {
                make_operator() {}
            };

            template <class T, char, class O>
            struct half_apply
            {
                T &&lhs;
            };

            template <class Lhs, class Op>
            half_apply<Lhs, '*', Op> operator*(Lhs &&lhs, make_operator<Op>)
            {
                return {std::forward<Lhs>(lhs)};
            }

            template <class Lhs, class Op, class Rhs>
            auto operator*(half_apply<Lhs, '*', Op> &&lhs, Rhs &&rhs)
                -> decltype(named_invoke(std::forward<Lhs>(lhs.lhs), Op{}, std::forward<Rhs>(rhs)))
            {
                return named_invoke(std::forward<Lhs>(lhs.lhs), Op{}, std::forward<Rhs>(rhs));
            }
        }
        namespace in_helper
        {
            struct in_t : notstd::named_operator::make_operator<in_t>
            {
            };
            template <class T, class C>
            bool named_invoke(T &&t, in_t, C &&c)
            {
                return notstd::find(std::forward<C>(c), std::forward<T>(t));
            }
        }
        // in_helper::in_t in;
        extern const in_helper::in_t in; // Declare `in` here as extern, and define it in util.cc

        bool strfind(const std::string &substring, const std::string &str);
    }

    namespace vector
    {
        std::vector<std::vector<int>> splitVectorByDelimiter(const std::vector<int> &input, int delimiter);

    }

}

using util::notstd::in;
using util::py::print;

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