#Class takes in a configuration file and returns a python dictionary format.
#Also has checks to make sure types aren't dumb
import json

class ConfigParser(object):
    def __init__(self, configdir):
        self.configdir = configdir
        print("CONFIGDIRECTORY: " + str(self.configdir))

    def Load_JsonConfig(self,filename):
        '''Load a json file with the filename in current self.configdir'''
        with open("%s/%s"%(self.configdir,filename),"r") as f:
            configdict = json.load(f)
        return configdict
    
    def SaveConfiguration(self,config_dict, savedir,config_outname):
        saveconfigloc = savedir+"/"+config_outname
        with open(saveconfigloc,"w") as f:
            json.dump(config_dict, f, sort_keys=True,indent=4)
