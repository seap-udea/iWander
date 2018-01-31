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
# Run a given task saving log information
#########################################################
#!/bin/bash
<<COMMENT
Usage:

   bash bin/run.sh wanderer.exe

   Where: 
   wanderer.exe is the progam to run

   Output: standard output will be redirected to
   log/<progname>-<wanderer>.log and
   log/<progname>-<wanderer>-detailed.log

COMMENT

############################################################
#PARAMETERS
############################################################
WANDERER=$(grep "char Wanderer" iwander.conf |head -n 1 |awk -F"\"" '{print $2}')

############################################################
#INPUT PARAMETERS
############################################################
program=$1;shift
progname=$(echo $program |cut -f 1 -d'.')
if [ ! -e $progname.cpp ];then
    echo "The program '$progname.cpp' does not exist."
    exit 1
fi

eval $@

############################################################
#COMPILING
############################################################
echo "Compiling $program..."
rm $program
make $program
if [ $? -gt 1 ];then
    echo "There was an error compiling $program.  Please check the code."
    exit 1
fi

############################################################
#RUNNING
############################################################
logfile="log/$progname-$WANDERER.log"
logdet="log/$progname-$WANDERER-detailed.log"
echo "Running (the output will appear in a moment and it is redirected to '$logfile')..."
./$program 2> $logdet > >(tee $logfile)
