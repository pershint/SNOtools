//This will have the main functions for parsing the config file we will write
//Will be associated with functions defined in the header file
//
//
#include "ConfigParser.hh"
#include "Exceptions.hh"
#include <exception>
#include <stdexcept>

namespace configuration {

    //Methods for Cofiguration Parser class
  CoParser::CoParser(const std::string &fName)
  {
    this->fName = fName;
    ExtractKeys();
  }

  void CoParser::ExtractKeys()
  {
    std::ifstream file;
    file.open(fName.c_str());
    if (!file) {
      std::cout << "CFG: file " << fName << " cannot be found...\n" << std::endl;
      throw(std::exception());
    }
    std::string line;
    size_t lineNo=0;
    while (std::getline(file,line))
    {
      lineNo++;
      std::string temp = line;
      if (temp.empty())
          continue;
      removeComment(temp);
      if (onlyWhiteSpace(temp))
          continue;
      parseLine(temp,lineNo);
    }
  }

  bool CoParser::onlyWhiteSpace(const std::string &line) const
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
    converter convert;
    if (line.find('=') == line.npos){
      std::cout << "CFG: Couldn't find the separator on line: " <<
        convert.T_to_string(lineNo) << "\n" << std::endl;
      throw(the_ex);
    }
    if (!(validLine(line))){
      std::cout << "CFG: Bad format for line: " << convert.T_to_string(lineNo) 
          << "\n" << std::endl;
      throw(the_ex);
    }
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
    else{
      std::cout << "CFG: Adding a key that already exists!  Aborting...\n"<<std::endl;
      throw(the_ex);
    }
  }

  void CoParser::extractKey(std::string &key, size_t const &sepPos, const std::string &line) const
  {
    key = line.substr(0, sepPos);
    if (key.find('\t') != line.npos || key.find(' ') != line.npos)
     key.erase(key.find_first_of("\t "));
  }

  void CoParser::extractValue(std::string &value, size_t const &sepPos, const std::string &line) const
  {
    value = line.substr(sepPos+1);
    value.erase(0,value.find_first_not_of("\t "));
    value.erase(value.find_last_not_of("\t ") + 1);
  }

  bool CoParser::keyExists(const std::string &key) const
  {
    return contents.find(key) != contents.end();
  }

} //Namespace configuration
