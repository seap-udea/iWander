#########################################################
#    _ _       __                __         		#
#   (_) |     / /___ _____  ____/ /__  _____		#
#  / /| | /| / / __ `/ __ \/ __  / _ \/ ___/		#
# / / | |/ |/ / /_/ / / / / /_/ /  __/ /    		#
#/_/  |__/|__/\__,_/_/ /_/\__,_/\___/_/     		#
# Dynamics of Interestellar Wanderers			#
# Jorge I. Zuluaga et al. [)] 2017			#
# http://github.com/seap-udea/iWander.git		#
#########################################################
# Generate the progenitor candidates CANDIDATES.md
#########################################################
#!/usr/bin/env python
############################################################
#PACKAGES
############################################################
from iwander import *

############################################################
#READ GENERAL CONFIGURATION FILE
############################################################
conf=readConf("iwander.conf");

############################################################
#PAST OR FUTURE
############################################################
qpast=1
if conf["duration"][0]=="-":
    candtxt="Progenitor candidates"
else:
    qpast=0
    candtxt="Future close encounters"

if qpast:
    pcols="Ppos|Pvmed|Pdist|Pprob|"
    lcols="--|--|--|--|"
    suffix="Past"
else:
    pcols="Ppos|"
    lcols="--|"
    suffix="Future"

############################################################
#LOAD PROGENITORS
############################################################
try:
    progenitors=pd.read_csv("scratch/progenitors-%s.csv"%conf["WANDERER"])
except FileNotFoundError:
    print("No analysis has been found for wanderer %s"%conf["WANDERER"])
    exit(1)

print("Number of progenitors:",len(progenitors))

#SELECT COLUMNS
sorting=conf["Sorting"].split(",")

progsort=progenitors.sort_values(by=sorting[0],ascending=eval(sorting[1]))

# MD Table
f=open("scratch/CANDIDATES-%s-%s.md"%(conf["WANDERER"],suffix),"w")
f.write("""# %s of %s

[![arXiv](http://img.shields.io/badge/arXiv-1711.09397-orange.svg?style=flat)](http://arxiv.org/abs/1711.09397)

_Latest update_: ``%s``

|HIP/TYCHO|Name|tmin|dmin|vrel|tmin|dmin|vrel|%s
|--|--|--|--|--|--|--|--|%s
"""%(candtxt,conf["WANDERER_NAME"],time.strftime("%c"),pcols,lcols))

i=1
for index in progsort.index:
    p=progenitors.loc[index]
    row=""
    if str(p.hip)=='nan':
        if str(p.tycho2_id)=='nan':sid='--'
        else:sid="TYC "+str(p.tycho2_id)
    else:sid="HIP "+str(int(p.hip))
    
    bf="";mbf=""

    if p.Psurmed!=0:
        if np.log10(p.Psurmed)<-20:continue
    else:
        continue
    
    simbad=str(p.name_simbad).replace('nan','--').replace('_',' ')
    if simbad!="--":
        simbad_ns=simbad.replace(" ","%20")
        simbad="[%s](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s)"%(simbad,simbad_ns)
    
    row+=r"| %s %s | %s %s | "%(bf,sid,
                              bf,simbad,
                            )
    row+=r"%s%.1f | %s%.1f | %s%.0f | "%(mbf,p.nomtmin/1e6,mbf,p.nomdmin,mbf,p.nomvrel)
    row+=r"%s[%.1f,%.1f,%.1f] | %s[%.1f,%.1f,%.1f] | %s[%.0f,%.0f,%.0f] | "%(mbf,p.tminl/1e6,p.tminmed/1e6,p.tminu/1e6,
                                                                             mbf,p.dminl,p.dminmed,p.dminu,
                                                                             mbf,p.vrell,p.vrelmed,p.vrelu)
    logPsurmed="%s%.1f"%(mbf,np.log10(p.Psurmed)) if p.Psurmed>0 else '--'
    if qpast:
        logPvelmed="%s%.1f"%(mbf,np.log10(p.Pvelmed)) if p.Pvelmed>0 else '--'
        logPdist="%s%.1f"%(mbf,np.log10(p.Pdist)) if p.Pdist>0 else '--'
        logPprob="%s%.1f"%(mbf,np.log10(p.Pprob)) if p.Pprob>0 else '--'
        row+=r"%s | %s | %s | %s |"%(logPsurmed,logPvelmed,logPdist,logPprob)
    else:
        row+=r"%s |"%(logPsurmed)
    f.write(row+"\n")
    i+=1
