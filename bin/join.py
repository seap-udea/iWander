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
#READ PROGENITOR FILES
############################################################
progenitors=pd.DataFrame()
for progfile in argv[1:]:
    print("\tMerging '%s'..."%progfile)
    try:
        data=pd.read_csv(progfile)
        progenitors=progenitors.append(data,ignore_index=True)
    except:
        print("\t\tNo data")

progfile="scratch/progenitors-%s.csv"%conf["Wanderer"]
print("Saving joined progenitors to '%s'..."%progfile)
progenitors.to_csv(progfile,index=False)
print("Number of progenitors:",len(progenitors))
