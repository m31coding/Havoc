#include <map>
#include <fstream>
#include <sstream>
#include "read_constants.h"

using namespace std;

// initialize key values
map <string, string> ReadConstants::keyValues = map<string, string>();

// read constants from file
void ReadConstants::readConstants(const char* filename)
{
    ifstream stream(filename);

    if (!stream.good())
    {
        fprintf(stderr, "ERROR in ReadConstants::readConstants(): can't read file\n");
        exit(1);
    }
    string line;

    while (getline(stream, line))
    {
        parseLine(line);
    }
}

// parse a line
void ReadConstants::parseLine(string& line)
{
    char key[256];
    char value[256];
    char* current = NULL;
    int currentIndex = 0;
    bool newToken = true;

    for (string::size_type i = 0; i < line.size(); ++i)
    {
        char c = line[i];

        if (c == ' ' || c == '\t')
        {
            newToken = true;
            continue;
        }

        if (c == '#')
        {
            return;
        }

        if (newToken)
        {
            if (current == NULL)
            {
                current = key;
                current[currentIndex++] = c;
            }
            else if (current == key)
            {
                current[currentIndex++] = '\0';
                current = value;
                currentIndex = 0;
                current[currentIndex++] = c;
            }
            else
            {
                break;
            }
        }
        else
        {
            current[currentIndex++] = c;
        }

        newToken = false;
    }

    if (current == value)
    {
        current[currentIndex++] = '\0';
        keyValues[key] = value;
    }
}

// check if key exists
bool ReadConstants::keyExists(const string& key)
{
    return keyValues.count(key) == 1;
}

// get a string value
string ReadConstants::GetString(const string& key)
{
    checkKey(key);
    return keyValues[key];
}

// get an integer value
int ReadConstants::GetInt(const string& key)
{
    checkKey(key);
    return strtol(keyValues[key].c_str(), NULL, 0);
}

// get an unsigned integer value
unsigned int ReadConstants::GetUnsignedInt(const string& key)
{
    checkKey(key);
    return strtoul(keyValues[key].c_str(), NULL, 0);
}

// get a double value
double ReadConstants::GetDouble(const string& key)
{
    checkKey(key);
    return strtod(keyValues[key].c_str(), NULL);
}

// get a boolean value
bool ReadConstants::GetBool(const string& key)
{
    checkKey(key);

    if (key == "true" || key == "on" || key == "yes")
    {
        return true;
    }

    if (key == "false" || key == "off" || key == "no")
    {
        return false;
    }

    return strtol(keyValues[key].c_str(), NULL, 0);
}

// print error and exit if the key does not exist
void ReadConstants::checkKey(const string& key)
{
    if (!keyExists(key))
    {
        fprintf(stderr, "The key %s does not exist", key.c_str());
        exit(0);
    }
}

// free memory
void ReadConstants::del()
{
    keyValues.clear();
}