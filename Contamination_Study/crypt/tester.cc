#include "ConfigParser.hh"

using namespace configuration;

int main(int argc, char* argv[]){
  configuration::CoParser cfgp("../config/cuts_default.ini");
  if (cfgp.keyExists("cut_DCmask")){
    std::cout << "the pathological mask exists!\n";
    //Grabbing values from the config file now
    int pathval = cfgp.getValueOfKey<int>("cut_DCmask");
    std::cout << "pathval: " << pathval << std::endl;
  }
}
 
