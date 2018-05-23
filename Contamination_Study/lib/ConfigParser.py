#Class takes in a configuration file and returns a python dictionary format.
#Also has checks to make sure types aren't dumb
import json

class ConfigParser(object):
    def __init__(self, configdir):
        self.configdir = configdir
        print("CONFIGDIRECTORY: " + str(self.configdir))

    def Parse_ini(self,filename):
        '''Parse a .ini file format in current defined self.configidr'''
        configdict = {}
        with open(self.configdir,"r") as f:
            for line in f:
                if line.startswith(";"):
                    continue
                if line.startswith(" "):
                    continue
                if not line.strip(): #line only has whitespace
                    continue
                configline = line.lstrip(";")
                keyval = line.split("=")
                key, value = keyval[0], keyval[1]
                value.split(";",1)[0].rstrip(" ")
                configdict[str(key)]=value.rstrip('\n')
        return configdict
        
    def Load_JsonConfig(self,filename):
        '''Load a json file with the filename in current self.configdir'''
        with open("%s/%s"%(self.configdir,filename),"r") as f:
            configdict = json.load(f)
        return configdict
    
    def SaveConfiguration(self,config_dict, savedir,config_outname):
        saveconfigloc = savedir+"/"+config_outname
        with open(saveconfigloc,"w") as f:
            json.dump(config_dict, f, sort_keys=True,indent=4)
