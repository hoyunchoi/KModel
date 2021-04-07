#! /bin/bash

networkSize=$1
SE=$2   #0.11
XQX=$3  #0.99
T=$4
E=$5
machine=$6
coreNum=$7

spg run ${machine} ./SA_KM.sh ${networkSize} ${SE} ${XQX} ${coreNum} ${T} ${E} ${coreNum}