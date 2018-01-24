//This will have the main functions for parsing the config file we will write
//Will be associated with functions defined in the header file
//
//
#include "ConfigParser.hh"
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <typeinfo>

namespace configuration {

  class converter  {
    template <typename T>
    static std::string T_to_string(T const &val)
    {
      std::ostringstream ostr;
      ostr << val;
      return ostr.str();
    }

    template <typename T>
    static T string_to_T(std::string const &val)
    {
      std::istringstream istr(val);
      T returnVal;
      if (!(istr >> returnVal))
        exitWithError("CFG: Not a valid " + (std::string)typeid(T).name() + " received!\n");
      return returnVal;
    }
  };


  template <typename ValueType>
  ValueType CoParser::getValueOfKey(const std::string &key, ValueType const
          &defaultValue = ValueType()) const
  {
    if (!keyExists(key))
      return defaultValue;

    return Convert::string_to_T<ValueType>(contents.find(key)->second);
  }

  inline CoParser::CoParser(const std::string &fName)
  {
    this->fName = fName;
    ExtractKeys();
  }

  void CoParser::ExtractKeys()
  {
    std::ifstream file;
    file.open(fName.c_str());
    if (!file)
      exitWithError("CFG: file" + fName.c_str() + "cannot be found...\n");
    std::string line;
    size_t lineNo=0;
    while (std::getline(file,line))
    {
      lineNo++;
      std::string temp = line;
      if temp.empty()
          continue;
      removeComment(temp);
      if onlyWhiteSpace(temp)
          continue;
      parseLine(line);
    }
  }

  bool CoParser::onlyWhitespace(const std::string &line) const
  {
    return (line.find_first_not_of(' ') == line.npos);
  }


  void CoParser::removeComment(std::string &line) const
  {
    if (line.find(';') != line.npos)
      line.erase(line.find(';'));
  }

  bool CoParser::validLine(const std::string &line) const
  {
    std::string temp = line;
    temp.erase(0, temp.find_first_not_of("\t "));
    if (temp[0] == '=')
      return false;
    for (size_t i = temp.find('=') + 1; i < temp.length(); i++){
      if (temp[i] != ' ')
        return true;
    } 
    return false;
  }

  void CoParser::parseLine(const std::string &line, size_t const lineNo)
  {
    if (line.find('=') == line.npos){
      exitWithError("CFG: Couldn't find the separator on line: " +
          converter::T_to_string(lineNo) + "\n");
    }
    if !(validLine(line))
        exitWithError("CFG: Bad format for line: " + converter::T_to_string(lineNo) +"\n");
    extractContents(line);
  }

  void CoParser::extractContents(const std::string &line)
  {
    std::string temp = line;
    temp.erase(0,temp.find_first_not_of("\t "));
    size_t sepPos = temp.find('=');
    //Extract the key and value pair
    std::string key, value;
    extractKey(key, sepPos, temp);
    extractValue(value, sepPos, temp);

    if (!keyExists(key))
        contents.insert(std::pair<std::string, std::string>(key,value));
    else
        exitWithError("CFG: Adding a key that already exists!  Aborting...\n");
  }

  void CoParser::extractKey(std::string &key, size_t const &sepPos, const std::string &line) const
  {
    key = line.substr(0, sepPos);
    if (key.find('\t') != line.npos || key.find(' ') != line.npos)
     key.erase(key.find_first_of("\t "));
  }

  void CoParser::extractValue(std::string &value, size_t const &sepPos, const std::string &line) const
  {
    value = line.substr(sepPos);
    value.erase(value.find_first_not_of("\t "));
    value.erase(value.find_last_not_of("\t ") + 1);
  }

  bool CoParser::exitWithError(std::string &error)
  {
    std::cout << error << std::endl;
    std::cin.ignore();
    std::cin.get();

    exit(EXIT_FAILURE);
  }

  bool CoParser::keyExists(std::string &key)
  {
    return contents.find(key) != contents.end();
  }
} //Namespace configuration
