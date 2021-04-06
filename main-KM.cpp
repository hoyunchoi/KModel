#include <chrono>
#include <iostream>
#include <random>  //random_device()
#include <string>

#include "../library/stringFormat.hpp"

#include "KM.hpp"

const double getEnergy(const std::vector<std::vector<unsigned>>&);

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
        int coreNum = std::stoi(argv[11]);
    */
    //* Base directory data will be saved
    const std::string dataDirectory = "../data/KM/";

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

    //* Generate K-Model and path for data
    KM model(network, rates, randomEngine);
    const std::string networkDirectory = model.directoryName(networkType, networkSize, meanDegree);
    CSV::generateDirectory(dataDirectory + networkDirectory);
    const std::string fileName = model.fileName(randomEngineSeed);

    //* Run and save
    auto start = std::chrono::system_clock::now();
    //*----------------------------------------------------------------------------------
    const std::vector<std::vector<unsigned>> simulatedData = model.sync_run(deltaT, maxDate);
    if (simulatedData.size() < maxDate){
        std::cout << "Simulated finished before reacing current time.\n";
        CSV::write(dataDirectory + networkDirectory + fileName, simulatedData);
    }
    else{
        std::cout << "Successfully finished.\n";
        CSV::write(dataDirectory + networkDirectory + fileName, simulatedData);
        const double energy = getEnergy(simulatedData);
        std::ofstream energyFile("energyList.txt", std::ios_base::app);
        energyFile << networkDirectory << fileName << ": " << energy << "\n";
    }
    //*----------------------------------------------------------------------------------
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream timeLog("time.log", std::ios_base::app);
    timeLog << networkDirectory << fileName << ": " << sec.count() << " second \n";

    return 0;
}


const double getEnergy(const std::vector<std::vector<unsigned>>& t_simulatedData){
    double RMSE = 0.0;

    //* Read real data
    std::vector<unsigned> realConfirmed;
    CSV::read("realConfirmed_208.txt", realConfirmed);
    const unsigned maxDate = realConfirmed.size();

    //* Get confirmed case from simulated data
    std::vector<unsigned> simulatedConfirmed;
    simulatedConfirmed.reserve(t_simulatedData.size());
    for (const std::vector<unsigned>& simulatedData : t_simulatedData) {
        simulatedConfirmed.emplace_back(simulatedData[7] + simulatedData[8] + simulatedData[10]);
    }

    //* Get RMSE from difference between simulation and real data
    for (unsigned date = 0; date < maxDate; ++date) {
        RMSE += std::pow((double)realConfirmed[date] - (double)simulatedConfirmed[date], 2.0);
    }
    std::cout << RMSE << std::endl;
    return std::sqrt(RMSE / (double)maxDate);
}