#!/bin/bash

# The following code will change a cst in all data files

# read from command line images and data dir
inList=$1
cst=$2
new=$3

#loop through file names
allF=`ls *$inList*`
echo $cst $new

# change name
for i in $allF
do
    newName=`echo $i | sed "s/$cst/$new/g"`  # required double quotes to expande to vars
    echo $newName
    mv $i $newName
done
