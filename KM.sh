#! /bin/bash

# Network input
networkSize=$1
meanDegree=10

# Rate input
SE=$2       #0.11
E_AI=0.39
pA=0.36
IQI=0.33
AR=0.11
QICR=0.08
XQX=$3      #0.09
tau=14

# Random seed and core num input
coreNum=$4

name=N${networkSize}M${meanDegree}SE${SE}XQX${XQX}C${coreNum}

g++ -O3 -march=native -flto -std=c++17 -o ./bin/${name}.out main-KM.cpp
./bin/${name}.out ${networkSize} ${meanDegree} ${SE} ${E_AI} ${pA} ${IQI} ${AR} ${QICR} ${XQX} ${tau} ${coreNum}
rm bin/${name}.out


