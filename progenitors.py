#!/usr/bin/env python
"""
Generate the progenitor candidates CANDIDATES.md
"""
############################################################
#PACKAGES
############################################################
from iwander import *

############################################################
#READ GENERAL CONFIGURATION FILE
############################################################
conf=readConf("iwander.conf");

############################################################
#LOAD PROGENITORS
############################################################
progenitors=pd.read_csv("progenitors-%s.csv"%conf["WANDERER"])
print("Number of progenitors:",len(progenitors))

#SELECT COLUMNS
progsort=progenitors.sort_values(by='Psurmed',ascending=False)
#progsort=progenitors.sort_values(by='nomdmin',ascending=True)
progsort[progsort.Pprob>0][["hip","tycho2_id","name_simbad","source",'Pprob', 'Psurmed', 'Pvelmed', 'Pdist', 'nomtmin', 'nomdmin',
       'nomvrel', 'mintmin', 'maxtmin', 'mindmin', 'maxdmin', 'minvrel',
       'maxvrel']]

# MD Table
f=open("CANDIDATES-%s.md"%conf["WANDERER"],"w")
f.write("""# Progenitor Candidates of 1I/2017 U1

[![arXiv](http://img.shields.io/badge/arXiv-1711.09397-orange.svg?style=flat)](http://arxiv.org/abs/1711.09397)

_Latest update_: ``%s``

|HIP/TYCHO|Name|tmin|dmin|vrel|tmin|dmin|vrel|Ppos|Pvmed|Pdist|Pprob|
|--|--|--|--|--|--|--|--|--|--|--|--|
"""%(time.strftime("%c")))

i=1
for index in progsort.index:
    p=progenitors.loc[index]
    row=""
    if str(p.hip)=='nan':
        if str(p.tycho2_id)=='nan':sid='--'
        else:sid="TYC "+str(p.tycho2_id)
    else:sid="HIP "+str(int(p.hip))
    
    bf="";mbf=""
    if np.log10(p.Psurmed)<-10:continue
    
    simbad=str(p.name_simbad).replace('nan','--').replace('_',' ')
    if simbad!="--":
        simbad_ns=simbad.replace(" ","%20")
        simbad="[%s](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s)"%(simbad,simbad_ns)
    
    row+=r"| %s %s | %s %s | "%(bf,sid,
                              bf,simbad,
                            )
    row+=r"%s%.1f | %s%.1f | %s%.0f | "%(mbf,p.nomtmin/1e6,mbf,p.nomdmin,mbf,p.nomvrel)
    row+=r"%s[%.1f,%.1f] | %s[%.1f,%.1f] | %s[%.0f,%.0f] | "%(mbf,p.mintmin/1e6,p.maxtmin/1e6,
                                                                    mbf,p.mindmin,p.maxdmin,
                                                                    mbf,p.minvrel,p.maxvrel)
    logPsurmed="%s%.1f"%(mbf,np.log10(p.Psurmed)) if p.Psurmed>0 else '--'
    logPvelmed="%s%.1f"%(mbf,np.log10(p.Pvelmed)) if p.Pvelmed>0 else '--'
    logPdist="%s%.1f"%(mbf,np.log10(p.Pdist)) if p.Pdist>0 else '--'
    logPprob="%s%.1f"%(mbf,np.log10(p.Pprob)) if p.Pprob>0 else '--'
    row+=r"%s | %s | %s | %s |"%(logPsurmed,logPvelmed,logPdist,logPprob)
    f.write(row+"\n")
    i+=1
