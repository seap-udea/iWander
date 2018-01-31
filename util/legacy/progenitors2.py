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

|#|HIP/TYCHO|Name|tmin|dmin|vrel|dmin|%s
|--|--|--|--|--|--|--|%s
"""%(candtxt,conf["WANDERER_NAME"],time.strftime("%c"),pcols,lcols))

i=1
for index in progsort.index:
    p=progsort.loc[index]

    if p.nomdmin>conf["dminMax"]:continue

    row=""
    if str(p.hip)=='nan':
        if str(p.tycho2_id)=='nan':sid='--'
        else:sid="TYC "+str(p.tycho2_id)
    else:sid="HIP "+str(int(p.hip))
    
    bf="";mbf=""

    simbad=str(p.name_simbad).replace('nan','--').replace('_',' ')
    if simbad!="--":
        simbad_ns=simbad.replace(" ","%20")
        simbad="[%s](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s)"%(simbad,simbad_ns)
    
    row+=r"| %d | %s %s | %s %s | "%(i,bf,sid,
                                bf,simbad)
    row+=r"%s%.3f | %s%.1f | %s%.0f | "%(mbf,p.nomtmin/1e6,mbf,p.nomdmin,mbf,p.nomvrel)
    row+=r"%s[%.2f,%.2f,%.2f] |"%(mbf,p.dminl,p.dminmed,p.dminu)
    Ppos="%s%.1f"%(mbf,p.Ppos)
    if qpast:
        Pvel="%s%.1f"%(mbf,p.Pvel)
        Pposvel="%s%.1f"%(mbf,p.Pposvel)
        Pdist="%s%.1f"%(mbf,p.Pdist)
        IOP="%s%.1f"%(mbf,p.IOP)
        row+=r"%s | %s | %s | %s | %s |"%(Ppos,Pvel,Pposvel,Pdist,IOP)
    else:
        row+=r"%s |"%(Ppos)
    print(row)
    f.write(row+"\n")
    i+=1

# Generate latex table
i=1
f=open("scratch/CANDIDATES-%s-%s.tex"%(conf["WANDERER"],suffix),"w")
f.write(r"""\begin{table*}
  \centering
  \scriptsize
  \begin{tabular}{ll|ccc|c|ccc}
  \hline
  \multicolumn{2}{c|}{ID} &  \multicolumn{3}{c|}{Nominal} & \multicolumn{1}{c|}{Range} & \multicolumn{3}{c}{$\log P$}  \\ \hline
  HIP/TYCHO & Other & $t\sub{min}$ & $d\sub{min}$ & $v\sub{rel}$ & $d\sub{min}$ & $P\sub{pos,vel}$ & $P\sub{dist}$ & IOP \\
            &       & Myr & pc & km/s & pc & & & \\
  \hline\hline
""")
for index in progsort.index:
    p=progenitors.loc[index]
    if p.nomdmin>conf["dminMax"]:continue

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
    
    row+=r"%s %s & %s %s & "%(bf,sid,
                              bf,simbad.replace('nan','--').replace('_',' '),
                            )
    row+=r"$%s{%.3f}$ & $%s{%.1f}$ & $%s{%.0f}$ & "%(mbf,p.nomtmin/1e6,mbf,p.nomdmin,mbf,p.nomvrel)
    row+=r"$%s{[%.2f,%.2f,%.2f]}$ & "%(mbf,p.dminl,p.dminmed,p.dminu)
    Pposvel="$%s{%.1f}$"%(mbf,p.Pposvel)
    Pdist="$%s{%.1f}$"%(mbf,p.Pdist)
    IOP="$%s{%.1f}$"%(mbf,p.IOP) 
    row+=r"%s & %s & %s \\"%(Pposvel,Pdist,IOP)
    f.write(row+"\n")
    print(row)
    i+=1
f.write(r"""\hline
  \end{tabular}
\caption{Interstellar origin probability (IOP) for a selected group of nearby stars.
\label{tab:Progenitors}}
\end{table*}
""")
