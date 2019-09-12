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
NSPLIT=0
MAXPROC=6
Wanderer=$(grep "char Wanderer" iwander.conf |head -n 1 |awk -F"\"" '{print $2}')
FILE="scratch/candidates-$Wanderer.csv"

############################################################
#SUBSTITUTE VARIABLES WITH COMMAND LINE OPTIONS
############################################################
eval $@

if [ ! -e $FILE ];then
    echo "You must run fist 'encounters.exe'"
    exit 1
fi

############################################################
#CHECK IF NSPLIT=1
############################################################
if [ $NSPLIT -eq 0 ];then 
    bash bin/run.sh probability.exe
    exit 0
fi

############################################################
#CREATE A CONFIGURATION FILE
############################################################
ndata=$(cat $FILE |wc -l)
nperproc=$((ndata/NSPLIT))
{
    echo -e "FILE=$FILE"
    echo -e "NDATA=$ndata"
    echo -e "Wanderer=$Wanderer"
    echo -e "MAXPROC=$MAXPROC"
    echo -e "NSPLIT=$NSPLIT"
} > scratch/probability-$Wanderer.conf

echo "Running parallel calculation of probabilities:"
cat scratch/probability-$Wanderer.conf
echo -e "OBJECTS PER PROC. = $nperproc"
echo -n "Press enter to continue or ctrl-c to cancel..."
read

############################################################
#SPLIT CANDIDATES
############################################################
python3 bin/split.py $FILE $NSPLIT

############################################################
#PREPARE
############################################################
make probability.exe

############################################################
#RUN
############################################################
n=0
ncand=0
while [ 1 ]
do
    i=1
    while [ $i -le $MAXPROC ]
    do
	in=$(printf "%05d" $n)
	fname="$FILE.$in"
	if [ -e $fname ];then
	    if [ ! -e "log/done-$Wanderer.$in" ];then
		ncand=$((ncand+nperproc))
		cmd="./probability.exe $FILE $n"
		$cmd &> log/probability.log.$in &
		echo "Calculating probability for candidates '$fname' (when completed candidates analysed $ncand) (PID = $!)..."
		((i++))
	    else
		if [ -e "log/start-$Wanderer.$in" ];then
		    echo "Waiting that the process finishes for '$fname'..."
		else
		    echo "Probabilities already computed for candidates '$fname'..."
		fi
	    fi
	else
	    echo "No more files available..."
	    exit 0
	fi
	((n++))
    done
    wait
done
