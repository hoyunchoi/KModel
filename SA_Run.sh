#! /bin/bash

networkSize=$1  # 160000
SE=$2           # 0.11
XQX=$3          # 0.09
T=$4            # 800
E=$5            # 80
machine=$6
coreNum=$7

spg run ${machine} ./SA_KM.sh ${networkSize} ${SE} ${XQX} ${T} ${E} ${coreNum}