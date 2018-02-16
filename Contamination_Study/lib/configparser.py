#Class takes in a configuration file and returns a python dictionary format.
#Also has checks to make sure types aren't dumb

class ConfigParser(object):
    def __init__(self, configfileloc):
        self.configfile = configfileloc
        self.config_dict = {}

    def Parse(self):
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


