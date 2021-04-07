#! /bin/bash

networkSize=$1
SE=$2   #0.11
XQX=$3  #0.99
machine=$4
coreNum=$5

# spg run ${machine} ./KM.sh ${networkSize} ${SE} ${XQX} ${coreNum}
spg run ${machine} ./KM_Rt.sh ${networkSize} ${SE} ${XQX} ${coreNum}