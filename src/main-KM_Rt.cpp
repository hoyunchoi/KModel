#include <chrono>
#include <iostream>
#include <random>  //random_device()
#include <string>

#include "Networks.hpp"
#include "stringFormat.hpp"
#include "linearAlgebra.hpp"

#include "KM_Rt.hpp"
#include "fileName.hpp"

int main(int argc, char* argv[]) {
    //* Base directory data will be saved
    using namespace COVID;
    const std::string dataDirectory = "KM_Rt_data/";

    //* Random variables
    const int coreNum = std::stoi(argv[11]);
    const int randomEngineSeed = coreNum == -1 ? (std::random_device())() : coreNum;
    pcg32 randomEngine(randomEngineSeed);

    //* Network variables
    const std::string networkType = "ER";
    const unsigned networkSize = std::stoul(argv[1]);
    const unsigned meanDegree = std::stoul(argv[2]);
    const unsigned long long linkSize = networkSize * meanDegree / 2;
    const Network<unsigned> network = ER::generate(networkSize, linkSize, randomEngine);
    // const double degreeExponent = 2.5;
    // const Network network = SF::generate(network, linkSize, degreeExponent, randomEngine);

    //* K-Model variables
    std::vector<double> rates(8, 0.0);
    rates[0] = std::stod(argv[3]);  // S_E
    rates[1] = std::stod(argv[4]);  // E_AI
    rates[2] = std::stod(argv[5]);  // pA
    rates[3] = std::stod(argv[6]);  // I_QI
    rates[4] = std::stod(argv[7]);  // A_R
    rates[5] = std::stod(argv[8]);  // QI_CR
    rates[6] = std::stod(argv[9]);  // X_QX
    rates[7] = std::stod(argv[10]); // tau
    const double deltaT = 1e-1;

    //* Read real data
    std::vector<unsigned> realConfirmed;
    CSV::read("realData/confirmed_209.txt", realConfirmed);
    // const unsigned maxDate = realConfirmed.size();
    const unsigned maxDate = 400;

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
        // std::ofstream energyFile("energyList.txt", std::ios_base::app);
        // energyFile << networkDirectory << rateFileName << ": " << model.getEnergy(realConfirmed) << "\n";
        model.save(dataDirectory + networkDirectory + rateFileName);
    }
    else{
        std::cout << "Simulated finished before reaching current time.\n";
        model.save(dataDirectory + networkDirectory + rateFileName);

    }
    //*----------------------------------------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream timeLog("time.log", std::ios_base::app);
    timeLog << "KM_Rt->" << networkDirectory << rateFileName << ": " << sec.count() << " second\n";

    return 0;
}
