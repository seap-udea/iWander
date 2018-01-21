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
#!/bin/bash

############################################################
#PARAMETERS
############################################################
NSPLIT=100
MAXPROC=6
WANDERER=$(grep "char WANDERER" iwander.conf |head -n 1 |awk -F"\"" '{print $2}')

############################################################
#SUBSTITUTE VARIABLES WITH COMMAND LINE OPTIONS
############################################################
eval $@

if [ ! -e scratch/candidates-$WANDERER.csv ];then
    echo "You must run fist 'encounters.exe'"
    exit 1
fi
ndata=$(cat scratch/candidates-$WANDERER.csv |wc -l)
nperproc=$((ndata/NSPLIT))
echo "Running parallel calculation of probabilities:"
echo -e "\tWANDERER = $WANDERER"
echo -e "\tMAXPROC = $MAXPROC"
echo -e "\tNSPLIT = $NSPLIT"
echo -e "\tOBJECTS PER PROC. = $nperproc"
echo
echo -n "Press enter to continue or ctrl-c to cancel..."
read

############################################################
#SPLIT CANDIDATES
############################################################
python3.5 bin/split.py $NSPLIT

############################################################
#PREPARE
############################################################
make probability.exe

############################################################
#RUN
############################################################
n=1
while [ 1 ]
do
    i=0
    while [ $i -le $MAXPROC ]
    do
	in=$(printf "%05d" $n)
	fname="scratch/candidates-$WANDERER.csv.$in"
	if [ -e $fname ];then
	    if [ ! -e "log/done-$WANDERER.$in" ];then
		echo "Calculating probability for candidates '$fname'..."
		./probability.exe $n >> log/probability.log.$in & 
		((i++))
	    else
		echo "Probabilities already computed for candidates '$fname'..."
	    fi
	else
	    echo "No more file available..."
	    exit 0
	fi
	((n++))
    done
    wait
done
