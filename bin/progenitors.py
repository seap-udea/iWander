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
if conf["Duration"][0]=="-":
    candtxt="Progenitor candidates"
else:
    qpast=0
    candtxt="Future close encounters"

if qpast:
    pcols="Ppos|Pvel|Pposvel|Pdist|IOP|"
    lcols="--|--|--|--|--|"
    suffix="Past"
else:
    pcols="Ppos|"
    lcols="--|"
    suffix="Future"

############################################################
#LOAD PROGENITORS
############################################################
try:
    progenitors=pd.read_csv("scratch/progenitors-%s.csv"%conf["Wanderer"])
except FileNotFoundError:
    print("No analysis has been found for wanderer %s"%conf["Wanderer"])
    exit(1)

print("Number of progenitors:",len(progenitors))

#SELECT COLUMNS
sorting=conf["Sorting"].split(",")

progsort=progenitors.sort_values(by=sorting[0],ascending=eval(sorting[1]))

# MD Table
f=open("scratch/CANDIDATES-%s-%s.md"%(conf["Wanderer"],suffix),"w")
f.write("""# %s of %s

[![arXiv](http://img.shields.io/badge/arXiv-1711.09397-orange.svg?style=flat)](http://arxiv.org/abs/1711.09397)

_Latest update_: ``%s``

|#|Name|d(pc)|q|dmin(pc)|tmin(Myr)|vrel(km/s)|%s
|--|--|--|--|--|--|--|%s
"""%(candtxt,conf["Wanderer_Name"],time.strftime("%c"),pcols,lcols))

i=1
for index in progsort.index:
    p=progsort.loc[index]
    
    if p.nomdmin>conf["DminMax"]:continue
    if p.qastro<1:continue

    row=""
    if str(p.hip)=='nan':
        if str(p.tycho2_id)=='nan':sid='--'
        else:sid="TYC "+str(p.tycho2_id)
    else:sid="HIP "+str(int(p.hip))
    
    bf="";mbf=""

    simbad=str(p.name_simbad).replace('nan','--').replace('_',' ')
    if simbad!="--":
        simbad_ns=simbad.replace(" ","%20")
        simbad="([%s](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s))"%(simbad,simbad_ns)
    else:simbad=""

    d=AU/np.tan(p.parallax/(60*60*1000.0)*DEG)/PARSEC;
    row+=r"| %d | %s %s %s | %.1f | %d | "%(i,bf,sid,simbad,d,p.qastro)
    row+=r"%s%.2f [%.2f,%.2f,%.2f] |%s %.2f [%.2f,%.2f,%.2f] |%s %.1f [%.1f,%.1f,%.1f] | "%(mbf,p.nomdmin,p.dminl,p.dminmed,p.dminu,
                                                                                           mbf,p.nomtmin/1e6,p.tminl/1e6,p.tminmed/1e6,p.tminu/1e6,
                                                                                           mbf,p.nomvrel,p.vrell,p.vrelmed,p.vrelu)
    Ppos="%s%.1f"%(mbf,p.Ppos)
    if qpast:
        Pvel="%s%.1f"%(mbf,p.Pvel)
        Pposvel="%s%.1f"%(mbf,p.Pposvel)
        Pdist="%s%.1f"%(mbf,p.Pdist)
        IOP="%s%.1f"%(mbf,p.IOP)
        row+=r"%s | %s | %s | %s | %s |"%(Ppos,Pvel,Pposvel,Pdist,IOP)
    else:
        row+=r"%s |"%(Ppos)

    f.write(row+"\n")
    i+=1

# Generate latex table
i=1
f=open("scratch/CANDIDATES-%s-%s.tex"%(conf["Wanderer"],suffix),"w")
f.write(r"""\begin{table*}
\centering
\scriptsize
\begin{tabular}{llll|ccc|ccccc}
\hline
\multicolumn{4}{c|}{Basic properties} & 
\multicolumn{3}{c|}{Encounter conditions} & \multicolumn{5}{c}{$\log P$}  \\ \hline
\# & Name & $d_*$ (pc) & $q$ & 
$t\sub{min}$ (Myr)  & 
$d\sub{min}$ (pc)   & 
$v\sub{rel}$ (km/s) & 
$P\sub{pos}$ & $P\sub{vel}$ & $P\sub{pos,vel}$ & $P\sub{dist}$ & IOP \\
  \hline\hline
""")
for index in progsort.index:
    p=progenitors.loc[index]

    if p.nomdmin>conf["DminMax"]:continue
    if p.qastro<1:continue

    row=""
    if str(p.hip)=='nan':
        if str(p.tycho2_id)=='nan':sid='--'
        else:sid="TYC "+str(p.tycho2_id)
    else:sid="HIP "+str(int(p.hip))
    bf="";mbf=""

    simbad=str(p.name_simbad).replace('nan','--').replace('_',' ')
    if simbad!="--":
        simbad_ns=simbad.replace(" ","%20")
        simbad=r"\href{http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s}{%s}"%(simbad_ns,simbad)

    if simbad=="nan":simbad=""
    else:simbad="(%s)"%simbad.replace('_',' ')
    if simbad=="(--)":simbad=""

    d=AU/np.tan(p.parallax/(60*60*1000.0)*DEG)/PARSEC;
    Ppos="$%s{%.1f}$"%(mbf,p.Ppos)
    Pvel="$%s{%.1f}$"%(mbf,p.Pvel)
    Pposvel="$%s{%.1f}$"%(mbf,p.Pposvel)
    Pdist="$%s{%.1f}$"%(mbf,p.Pdist)
    IOP="$%s{%.1f}$"%(mbf,p.IOP) 

    row+=r"%s %d & %s%s & %.1f & %d & "%(bf,i,bf,sid,p.d,p.qastro)
    row+=r"$%s{%.2f}$ &"%(mbf,p.nomtmin/1e6)
    row+=r"$%s{%.2f}$ &"%(mbf,p.nomdmin)
    row+=r"$%s{%.1f}$ &"%(mbf,p.nomvrel)
    row+=r"%s & %s & %s & %s & %s \\"%(Ppos,Pvel,Pposvel,Pdist,IOP)
    f.write(row+"\n")
    
    row=""
    row+=r" & %s%s & & & "%(bf,simbad)
    row+=r"\tiny $%s[%.2f,%.2f,%.2f]$ &"%(mbf,p.tminl/1e6,p.tminmed/1e6,p.tminu/1e6)
    row+=r"\tiny $%s[%.2f,%.2f,%.2f]$ &"%(mbf,p.dminl,p.dminmed,p.dminu)
    row+=r"\tiny $%s[%.1f,%.1f,%.1f]$ &"%(mbf,p.vrell,p.vrelmed,p.vrelu)
    row+=r" & & & & \\"

    f.write(row+"\n\n")

    i+=1
f.write(r"""\hline
  \end{tabular}
\caption{Interstellar origin probability (IOP) for a selected group of nearby stars.
\label{tab:Progenitors}}
\end{table*}
""")
