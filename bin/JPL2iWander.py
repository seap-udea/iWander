#!/usr/bin/env
"""
   Convert from the JPL format to the configuration format of iWander
"""
copy="""
   From the JPL Small-Body Database Browser copy and paste starting in
   "Orbital Elements..." and finishing in the last line of the
   covariance matrix
"""
from sys import argv,exit
import re,os
from collections import OrderedDict

try:
    short=argv[1]
    body=argv[2]
    jplfile="%s.jpl"%short
except:
    print("You must provide a short and a long object name:")
    print("\tOumuamua \"1I/2017 U1 ('Oumuamua)\"")
    exit(1)

if not os.path.isfile(jplfile):
    print("You must first prepare the JPL file '%s' for object:"%jplfile)
    print(copy)
    exit(2)

coms="\tFor updated values see: https://ssd.jpl.nasa.gov/sbdb.cgi\n\n"
qcom=False
props=OrderedDict([["e","e"],["a","a"],["q","q"],["i","i"],
                   ["node","W"],["peri","w"],["M","M"],["n","n"],["tp","tp"]])
nom=OrderedDict()
err=OrderedDict()
uni=OrderedDict()
cov=OrderedDict()
print("Parsing JPL solution:")
for line in open(jplfile).readlines():
    line=line.strip()

    if "Orbital Elements" in line:
        m=re.search("Epoch\s*(\d+\.\d+)\s*\(([^)]+)\)\s+TDB",line)
        ini_to_jed=m.group(1)
        ini_to_date=m.group(2)
        print("\tinit_to_date = ",ini_to_date)
        print("\tinit_to_jed = ",ini_to_jed)
    if "Reference" in line:
        m=re.search("Reference:\s*(JPL\s*\d+)",line)
        reference_solution=m.group(1)
        print("\treference_solution = ",reference_solution)

    for prop in props.keys():
        if not re.search("^%s\s+"%prop,line) is None:
            parts=line.split()
            if len(parts)<6:
                #Read nominal and errors
                print("\tReading properties for %s:"%prop)
                nom[prop]=parts[1]
                if prop!="tp": 
                    err[prop]=parts[2]
                    try:uni[prop]=parts[3]
                    except:uni[prop]="adim."
                else:
                    err[prop]=""
                    uni[prop]="jed"
                print("\t\tNominal = ",nom[prop])
                print("\t\tError = ",err[prop])
                print("\t\tUnits = ",uni[prop])
            else:
                try:
                    value=float(parts[1])
                    #Read matrix
                    print("\tReading covariance matrix for %s:"%prop)
                    cov[prop]=parts[1:]
                    print("\t\tCovariance components = ",cov[prop][:3],"...")
                except:pass

    if "JED" in line:
        m=re.search("^\(([^)]+)\)\s*(\d+\.\d+)\s*JED",line)
        ini_tp_date=m.group(1)
        ini_tp_jed=m.group(2)
        print("\t\tPerihelion passtime = ",ini_tp_date)
        print("\t\tError in passtime = ",ini_tp_jed)

    if "obs. used" in line:qcom=True
    if "Additional" in line:qcom=False
    if qcom:
        if not re.search("\w+",line) is None:
            coms+="\t"+" ".join(line.split())+"\n"

print("\tComments:\n",coms)

#Generate configuration lines
lines="""/*
%s

%s
*/

//REFERENCE SOLUTION
SpiceChar ini_to_date[]="%s";
double ini_to_jed=%s;
SpiceChar reference_solution[]="%s";

//ORBITAL ELEMENTS WITH ERRORS
"""%(body,coms,ini_to_date,ini_to_jed,reference_solution)

#Orbital elements
for prop in props.keys():
    if prop=="tp":continue
    lines+="double ini_%s=%s,ini_d%s=%s;//%s\n"%(props[prop],
                                                 nom[prop],
                                                 props[prop],
                                                 err[prop],
                                                 uni[prop]
                                             )
#tp
lines+="""double ini_tp=%s;//%s
SpiceChar ini_tp_date[]=\"%s\";
double ini_tp_jed=%s;//%s"""%(nom["tp"],uni["tp"],ini_tp_date,
                              ini_tp_jed,uni["tp"])


#Covariance matrix
lines+="""

//COVARIANCE MATRIX
//ORDER IS: e,q,tp,W,w,i
double ini_cov[][6]={
"""
for prop in "e","q","tp","node","peri","i":
    covline="\t{"
    for value in cov[prop]:
        covline+=value+","
    covline=covline.strip(",")
    covline+="},\n"
    lines+=covline
lines=lines.strip(",")
lines+="};"

fconf="%s.conf"%short
print("Configuration file '%s' written"%fconf)
f=open(fconf,"w")
f.write(lines)
f.close()
print("Done.")

