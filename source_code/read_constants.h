#ifndef READCONSTANTS_H
#define READCONSTANTS_H

#include <map>

/**
a singleton class for reading constants from file
**/
class ReadConstants
{
public:

    /// read constants from file
    static void readConstants(const char* filename);

    /// check if key exists
    static bool keyExists(const std::string& key);

    /// get a string value
    static std::string GetString(const std::string& key);

    /// get an integer value
    static int GetInt(const std::string& key);

    /// get an unsigned integer value
    static unsigned int GetUnsignedInt(const std::string& key);

    /// get a double value
    static double GetDouble(const std::string& key);

    /// get a boolean value
    static bool GetBool(const std::string& key);

    /// free memory
    static void del();

private:

    /// key values map
    static std::map <std::string, std::string> keyValues;

    /// parse a line
    static void parseLine(std::string& line);

    /// print error and exit if the key does not exist
    static void checkKey(const std::string& key);
};

#endif