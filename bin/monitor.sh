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
# Monitor a parallel run
#########################################################
#!/bin/bash

############################################################
#PARAMETERS
############################################################
WANDERER=$(grep "char WANDERER" iwander.conf |head -n 1 |awk -F"\"" '{print $2}')
FILECONF="scratch/probability-$WANDERER.conf"
if [ ! -e scratch/probability-$WANDERER.conf ];then
    echo "No parallel job has been launched for this object"
    exit 1
fi
. scratch/probability-$WANDERER.conf

############################################################
#CHECK FOR FILES
############################################################
ncomp=$(ls log/done-$WANDERER.* 2>/dev/null |wc -l )
echo "Total number of candidates: $NDATA"
echo "Total number of tasks: $NSPLIT"
echo "Completed tasks: $ncomp"
echo "Stars completed per task:"

i=0
ntotal=0
for file in $(ls log/probability.log.*)
do
    in=$(echo $file |awk -F"." '{print $NF}')
    nstar=$(grep "Computing probabilities for candidate star" $file |tail -n 1 |awk '{print $6}')
    echo -e "\tNumber of stars completed for job $in: $nstar"
    ((i++))
    ntotal=$((ntotal+nstar))
done
echo "Total number of analysed stars: $ntotal"

