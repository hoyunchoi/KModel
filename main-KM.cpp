#include <chrono>
#include <filesystem>
#include <iostream>
#include <random>  //random_device()
#include <string>

#include "../library/stringFormat.hpp"
#include "KM.hpp"

int main(int argc, char* argv[]) {
    /*
        Input parameters
        unsigned networkSize = std::stoul(argv[1]);
        unsigned meanDegree = std::stoul(argv[2]);
        double S_E = std::stod(argv[3]);
        double E_A = std::stod(argv[4]);
        double E_I = std::stod(argv[5]);
        double I_QI = std::stod(argv[6]);
        double A_R = std::stod(argv[7]);
        double QI_CR = std::stod(argv[8]);
        double X_QX = std::stod(argv[9]);
        double tau = std::stod(argv[10]);
        int randomEngineSeed = std::stoi(argv[11]);
    */
    //* Random variables
    KM::randomEngineSeed = std::stoi(argv[11]);
    const int seed = KM::randomEngineSeed == -1 ? (std::random_device())() : KM::randomEngineSeed;
    pcg32 randomEngine(seed);

    //* Network variables
    const std::string networkType = "ER";
    const unsigned networkSize = std::stoul(argv[1]);
    const unsigned meanDegree = std::stoul(argv[2]);
    const unsigned linkSize = networkSize * meanDegree / 2;
    const Network network = ER::generate(networkSize, linkSize, randomEngine);
    // const double degreeExponent = 2.5;
    // const Network network = SF::generate(network, linkSize, degreeExponent, randomEngine);

    //* Rate variables
    KM::seedSize = 1;
    KM::S_E = std::stod(argv[3]);
    KM::E_AI = std::stod(argv[4]);
    KM::pA = std::stod(argv[5]);
    KM::E_A = KM::E_AI * KM::pA;
    KM::E_I = KM::E_AI * (1.0-KM::pA);
    KM::I_QI = std::stod(argv[6]);
    KM::A_R = std::stod(argv[7]);
    KM::QI_CR = std::stod(argv[8]);
    KM::X_QX = std::stod(argv[9]);
    KM::tau = std::stod(argv[10]);

    //* Model variables
    const double deltaT = 1e-1;
    KM::deletion = false;
    KM::directory = "../data/KM/" + networkType + "_N" + to_stringWithExponent((double)networkSize, 1) + ",M" + std::to_string(meanDegree) + "/";
    CSV::generateDirectory(KM::directory);

    auto start = std::chrono::system_clock::now();
    //* ---------------Initialize and Run------------------
    KM::initialize(network, randomEngine);
    // const double orderParameter = KM::asyncRun();
    const double orderParameter = KM::syncRun(deltaT);
    //*----------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    printf("order parameter: %0.2f, %0.10f second\n", orderParameter, sec.count());

    return 0;
}