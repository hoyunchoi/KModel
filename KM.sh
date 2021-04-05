#! /bin/bash

networkSize=20000
meanDegree=10
SE=$1       #0.11
E_AI=0.39
pA=0.36
IQI=0.33
AR=0.11
QICR=0.08
tau=14
XQX=$2      #0.09
seed=$3

name=N${networkSize}M${meanDegree}SE${SE}XQX${XQX}S${seed}

g++ -O3 -march=native -flto -std=c++17 -o bin/${name}.out main-KM.cpp
./bin/${name}.out ${networkSize} ${meanDegree} ${SE} ${E_AI} ${pA} ${IQI} ${AR} ${QICR} ${XQX} ${tau} ${seed}
rm bin/${name}.out


