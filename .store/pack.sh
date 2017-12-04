#!/bin/bash
basedir=$(pwd)
pack=$1
if [ "x$pack" == "x" ];then pack="pack";fi

storedir=".store"
confile="pack.conf"

if [ $pack = "pack" ];then
    echo "Packing..."
    rm -r $storedir/*--* &> /dev/null
    for file in $(cat $storedir/$confile |grep -v "#")
    do
	echo -e "\tfile $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	uname=$(echo $dname |sed -e s/\\//_/)
	sdir="$storedir/$uname--$fname"
	mkdir -p "$sdir"
	cd $sdir
	split -b 20000KB $basedir/$file $fname-
	cd - &> /dev/null
	git add "$storedir/$uname--$fname/*"
    done
    git add $storedir/*--*
else
    echo "Unpacking..."
    for file in $(cat $storedir/$confile |grep -v "#")
    do
	echo -en "\tUnpacking $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	uname=$(echo $dname |sed -e s/\\//_/)
	sdir="$storedir/$uname--$fname"
	if [ -e $dname/$fname ];then
	    echo "file already unpacked"
	else
	    cat "$sdir"/$fname-* > $dname/$fname
	    echo
	fi
    done
fi
