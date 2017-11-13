#!/bin/bash
pack=$1
if [ "x$pack" == "x" ];then pack="pack";fi

storedir=".store"
confile="pack.conf"

if [ $pack = "pack" ];then
    echo "Packing..."
    rm -r $storedir/*--*
    for file in $(cat $storedir/$confile)
    do
	echo -e "\tfile $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	sdir="$storedir/$dname--$fname"
	mkdir -p "$sdir"
	cd $sdir
	split -b 20000KB ../../$file $fname-
	cd - &> /dev/null
    done
    git add $storedir/*--*
else
    echo "Unpacking..."
    for file in $(cat $storedir/$confile)
    do
	echo -e "\tUnpacking $file..."
	fname=$(basename $file)
	dname=$(dirname $file)
	sdir="$storedir/$dname--$fname"
	cat "$sdir"/$fname-* > $dname/$fname
    done
fi
