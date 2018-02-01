import json

data = {}

with open("N16_ANALYSIS_final.txt","r") as f:
    for line in f:
        if "Run" in line:
            continue
        splitline = line.split(",")
        print(splitline)
        for j,entry in enumerate(splitline):
            if j==0:
                if len(splitline) < 10:
                    continue
                data[str(splitline[j].strip(" \t"))] = {}
                data[str(splitline[j].strip(" \t"))]["position"] = [float(splitline[j+1].strip(" \t")),float(splitline[j+2].strip(" \t")),float(splitline[j+3].strip(" \t"))]
print(data)
with open("tester.json","w") as f:
    json.dump(data,f,sort_keys=True,indent=4)
