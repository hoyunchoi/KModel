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
tau=14
XQX=$3     #0.09

# Simulated Annihilation input
T=$4
E=$5

# Core num input
coreNum=$6

name=N${networkSize}M${meanDegree}SE${SE}XQX${XQX}T${T}E${E}C${coreNum}

function debugBuild {
	g++ -std=c++17 -Wall -g -fsanitize=address\
        -I ${common} -I ${libDir}\
        -o ${binDir}/${name}\
	    ${srcDir}/main-SA_KM.cpp
}

function build {
	g++ -std=c++17 -O3 -flto -march=native\
        -I ${common} -I ${libDir}\
        -o ${binDir}/${name}\
		${srcDir}/main-SA_KM.cpp
}


#* Compile the source files
# build
debugBuild

#* Run
./bin/${name}.out ${networkSize} ${meanDegree} ${SE} ${E_AI} ${pA} ${IQI} ${AR} ${QICR} ${XQX} ${tau} ${T} ${E}
rm bin/${name}.out


