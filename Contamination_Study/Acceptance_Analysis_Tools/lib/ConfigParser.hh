//This will be the main header file for the ConfigParser class we will make
//Source: http://www.dreamincode.net/forums/topic/183191-create-a-simple-configuration-file-parser/
//
#include <string>
#include <map>

namespace configuration {

  class CoParser  {
    public:
      inline CoParser(const std::string &fName);
      bool keyExists(const std::string &key);
      template <typename ValueType>
      ValueType getValueOfKey(const std::string &key, ValueType const
          &defaultValue = ValueType());
    private:
      std::map<std::string, std::string> contents;
      std::string fName;

      //Main command; uses all others to fill in the contents map
      void ExtractKeys();
      void removeComment(std::string &line);  //Removes comments in lines
      bool onlyWhiteSpace(std::string &line); //True if line is blank space

      //Used to check each line is valid; if it is, contents extracted to contents
      void parseLine(const std::string &line, const size_t lineNo); //
      bool validLine(const std::string &line); //True if valid key=value format
      void extractContents(const std::string &line);
      void extractKey(std::string &key, size_t const &sepPos, const std::string &line);
      void extractValue(std::string &value, size_t const &sepPos, const std::string &line);
      //Error handling
      void exitWithError(std::string &error);
  }; 
};
