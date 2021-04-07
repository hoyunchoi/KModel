#include <chrono>
#include <iostream>
#include <random>  //random_device()
#include <string>

#include "../library/Networks.hpp"
#include "../library/stringFormat.hpp"
#include "../library/linearAlgebra.hpp"

#include "KM_Rt.hpp"
#include "fileName.hpp"

int main(int argc, char* argv[]) {
    //* Base directory data will be saved
    const std::string dataDirectory = "../data/KM_Rt/";

    //* Random variables
    const int coreNum = std::stoi(argv[11]);
    const int randomEngineSeed = coreNum == -1 ? (std::random_device())() : coreNum;
    pcg32 randomEngine(randomEngineSeed);

    //* Network variables
    const std::string networkType = "ER";
    const unsigned networkSize = std::stoul(argv[1]);
    const unsigned meanDegree = std::stoul(argv[2]);
    const unsigned linkSize = networkSize * meanDegree / 2;
    const Network network = ER::generate(networkSize, linkSize, randomEngine);
    // const double degreeExponent = 2.5;
    // const Network network = SF::generate(network, linkSize, degreeExponent, randomEngine);

    //* K-Model variables
    std::vector<double> rates(8, 0.0);  //* S_E, E_AI, pA, I_QI, A_R, QI_CR, X_QX, tau
    rates[0] = std::stod(argv[3]);
    rates[1] = std::stod(argv[4]);
    rates[2] = std::stod(argv[5]);
    rates[3] = std::stod(argv[6]);
    rates[4] = std::stod(argv[7]);
    rates[5] = std::stod(argv[8]);
    rates[6] = std::stod(argv[9]);
    rates[7] = std::stod(argv[10]);
    const double deltaT = 1e-1;
    const unsigned maxDate = 209;

    //* Generate K-Model with reproduction number and path for data
    KM_Rt model(network, rates, randomEngine);
    const std::string networkDirectory = networkName(networkType, networkSize, meanDegree);
    const std::string rateFileName = rateName(rates, randomEngineSeed);
    CSV::generateDirectory(dataDirectory + networkDirectory);

    //* Run and save
    auto start = std::chrono::system_clock::now();
    //*----------------------------------------------------------------------------------
    const bool finished = model.sync_run(deltaT, maxDate);
    if (finished){
        std::cout << "Successfully finished.\n";
    }
    else{
        std::cout << "Simulated finished before reaching current time.\n";
    }
    model.save(dataDirectory + networkDirectory + rateFileName);
    //*----------------------------------------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream timeLog("time.log", std::ios_base::app);
    timeLog << "KM_Rt->" << networkDirectory << rateFileName << ": " << sec.count() << " second\n";

    return 0;
}