#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "SA_KM.hpp"

int main(int argc, char* argv[]) {
    /*
        Input parameters
        unsigned networkSize = std::stoul(argv[1]);
        unsigned meanDegree = std::stoul(argv[2]);
        double S_E = std::stod(argv[3]);
        double E_AI = std::stod(argv[4]);
        double pA = std::stod(argv[5]);
        double I_QI = std::stod(argv[6]);
        double A_R = std::stod(argv[7]);
        double QI_CR = std::stod(argv[8]);
        double X_QX = std::stod(argv[9]);
        double tau = std::stod(argv[10]);
        const double initialTemperature = std::stod(argv[11]);
        const unsigned maxEnsemble = std::stoul(arvg[12]);
    */
    //* Base directory data will be saved
    const std::string dataDirectory = "../data/KM/";

    //* Network variables
    const unsigned networkSize = std::stoul(argv[1]);
    const unsigned meanDegree = std::stoul(argv[2]);
    // const double degreeExponent = 2.5;
    // const Network network = SF::generate(network, linkSize, degreeExponent, randomEngine);

    //* Parameters for Simulated Annihilation
    std::vector<double> initialRates(8, 0.0);  //* S_E, E_AI, pA, I_QI, A_R, QI_CR, X_QX, tau
    initialRates[0] = std::stod(argv[3]);
    initialRates[1] = std::stod(argv[4]);
    initialRates[2] = std::stod(argv[5]);
    initialRates[3] = std::stod(argv[6]);
    initialRates[4] = std::stod(argv[7]);
    initialRates[5] = std::stod(argv[8]);
    initialRates[6] = std::stod(argv[9]);
    initialRates[7] = std::stod(argv[10]);
    const double initialTemperature = std::stod(argv[11]);
    const unsigned maxEnsemble = std::stoul(argv[12]);

    //* Generate SA_KM model
    const double deltaT = 1e-1;
    SA_KM model(networkSize, meanDegree, initialRates, initialTemperature);

    //* Run and save
    auto start = std::chrono::system_clock::now();
    //*----------------------------------------------------------------------------------
    model.sync_run(maxEnsemble, deltaT);
    //*----------------------------------------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream timeLog("SA_time.log", std::ios_base::app);
    timeLog << "N" << networkSize << ",M" << meanDegree << ",E" << maxEnsemble << ": " << sec.count() << " second \n";

    return 0;
}