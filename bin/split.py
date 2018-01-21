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
#READ NUMBER OF PARTS ON 
############################################################
try:nparts=int(argv[1])
except:nparts=1

############################################################
#READ LIST OF CANDIDATES
############################################################
cname="scratch/candidates-%s.csv"%conf["WANDERER"]
data=pd.read_csv(cname)
ndata=len(data)
print("Total number of candidates:",ndata)
print("Total number of pars:",nparts)

############################################################
#SPLIT CANDIDATES
############################################################
dndata=int((1.*ndata)/nparts)
for i in range(nparts):
    ipart=i*dndata
    epart=ipart+dndata
    if i==(nparts-1):epart=ndata
    subdata=data.iloc[ipart:epart]
    spart=len(subdata)
    
    print("\tSaving Pack %d (%d,%d = %d)"%(i,ipart,epart,spart))
    subdata.to_csv(cname+".%05d"%i,index=False)

