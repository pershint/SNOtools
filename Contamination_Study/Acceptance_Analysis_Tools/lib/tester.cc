#include "ConfigParser.cc"
#include "ConfigParser.hh"


int main(int argc, char* argv[]){
    configuration::CoParser cfgp("../config/cuts.ini");
  if (cfgp.keyExists("path_DCmask")){
    std::cout << "the pathological mask exists!\n";
    //Grabbing values from the config file now
    int pathval = cfgp.getValueOfKey<int>("path_DCmask");
    cout << "pathval: " << pathval << endl;
  }
}
 
