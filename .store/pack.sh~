#!/bin/bash
pack=$1
if [ "x$pack" == "x" ];then pack="pack";fi

if [ $pack = "pack" ];then
    echo "Packing..."
    rm -r Store/*--*
    for file in $(cat Store/store.conf)
    do
	echo -e "\tfile $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	sdir="Store/$dname--$fname"
	mkdir -p "$sdir"
	cd $sdir
	split -b 20000KB ../../$file $fname-
	cd - &> /dev/null
    done
    git add Store/*--*
else
    echo "Unpacking..."
    for file in $(cat Store/store.conf)
    do
	echo -e "\tUnpacking $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	sdir="Store/$dname--$fname"
	cat "$sdir"/$fname-* > $dname/$fname
    done
fi
