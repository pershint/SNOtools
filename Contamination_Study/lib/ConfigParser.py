#Class takes in a configuration file and returns a python dictionary format.
#Also has checks to make sure types aren't dumb
import json

class ConfigParser(object):
    def __init__(self, configfileloc):
        self.configfile = configfileloc
        print("CONFIGFILE: " + str(self.configfile))
        self.config_dict = {}

    def Parse_ini(self):
        configdict = {}
        with open(self.configfile,"r") as f:
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
        self.config_dict = configdict
        return configdict
        
    def Load_JsonConfig(self):
        with open(self.configfile,"r") as f:
            configdict = json.load(f)
        self.config_dict = configdict
        return configdict

    def SaveConfiguration(self,config_dict, savedir,config_outname):
        saveconfigloc = savedir+"/"+config_outname
        with open(saveconfigloc,"w") as f:
            json.dump(config_dict, f, sort_keys=True,indent=4)
