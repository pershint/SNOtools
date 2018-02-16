//This will be the main header file for the ConfigParser class we will make
//Source: http://www.dreamincode.net/forums/topic/183191-create-a-simple-configuration-file-parser/
//
//LESSONS LEARNED:
// - Easiest way to deal with templates is to define explicitly
//   In the header file itself
// - Default variables can be defined in .hh or .cc, but not both
// - const definitions are a safeguard against writing to variables when no write
//   actions should be happening
// - A reminder that giving the address of a type in memory means you will read/write
//   to that addresse's contents.

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <typeinfo>
#include <cstdlib>

namespace configuration {

  class converter  {
    public:
      template <typename T>
        std::string T_to_string(T const &val)
        {
          std::ostringstream ostr;
          ostr << val;
          return ostr.str();
        }
        
        template <typename T>
        T string_to_T(std::string const &val)
        {
          std::istringstream istr(val);
          T returnVal;
          if (!(istr >> returnVal))
            exitWithError("CFG: Not a valid type id received!\n");
          return returnVal;
        }

      
        void exitWithError(std::basic_string<char> error)
        {
          std::cout << error << std::endl;
          std::cin.ignore();
          std::cin.get();
          exit(EXIT_FAILURE);
        }
  };

  class CoParser  {
    public:
      CoParser(const std::string &fName);
      bool keyExists(const std::string &key) const;
      template <typename ValueType>
      ValueType getValueOfKey(const std::string &key, const ValueType
              &defaultValue = ValueType()) const
      {
        if (!keyExists(key)){
          return defaultValue;
        }
        converter convert;
        return convert.string_to_T<ValueType>(contents.find(key)->second);
      }

    private:
      std::map<std::string, std::string> contents;
      std::string fName;

      //Main command; uses all others to fill in the contents map
      void ExtractKeys();
      void removeComment(std::string &line) const;  //Removes comments in lines
      bool onlyWhiteSpace(const std::string &line) const; //True if line is blank space

      //Used to check each line is valid; if it is, contents extracted to contents
      void parseLine(const std::string &line, const size_t lineNo); //
      bool validLine(const std::string &line) const; //True if valid key=value format
      void extractContents(const std::string &line);
      void extractKey(std::string &key, size_t const &sepPos, const std::string &line) const;
      void extractValue(std::string &value, size_t const &sepPos, const std::string &line) const;
      //Error handling
      void exitWithError(std::basic_string<char> error);
  }; 
};
