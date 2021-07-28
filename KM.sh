#! /bin/bash

srcDir=src
libDir=lib
binDir=bin
common=../library

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

function debugBuild {
	g++ -std=c++17 -Wall -g -fsanitize=address\
        -I ${common} -I ${libDir}\
        -o ${binDir}/${name}\
	    ${srcDir}/main-KM.cpp
}

function build {
	g++ -std=c++17 -O3 -flto -march=native\
        -I ${common} -I ${libDir}\
        -o ${binDir}/${name}\
		${srcDir}/main-KM.cpp
}

#* Compile the source files
# build
debugBuild

#* Run
./${binDir}/${name} ${networkSize} ${meanDegree} ${SE} ${E_AI} ${pA} ${IQI} ${AR} ${QICR} ${XQX} ${tau} ${coreNum}
rm ${binDir}/${name}
