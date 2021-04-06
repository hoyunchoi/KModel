#pragma once

#include <random>
#include <vector>
#include <string>
#include <cmath>

#include "../library/CSV.hpp"
#include "../library/Networks.hpp"
#include "../library/pcg_random.hpp"

#include "KM.hpp"

struct SA_KM {
    //* Member variables
   private:
    unsigned m_networkSize;
    unsigned m_meanDegree;
    std::vector<double> m_rates;  //* S_E, E_AI, pA, I_QI, A_R, QI_CR, X_QX, tau
    double m_temperature;
    pcg32 m_seedEngine;
    std::uniform_real_distribution<double> m_probDistribution;
    std::uniform_int_distribution<int> m_randomEngineSeedDistribution;
    std::vector<unsigned> m_realConfirmed;
    unsigned m_maxDate;
    double m_energy = 1e6;

    //* Member functions
   private:
    const double m_getEnergy(const std::vector<std::vector<unsigned>>&) const;
    const std::vector<double> m_perturbate(const double&) ;
   public:
    SA_KM() {}
    SA_KM(const unsigned&, const unsigned&, const std::vector<double>&, const double&);
    void sync_run(const unsigned&, const double&);
};

SA_KM::SA_KM(const unsigned& t_networkSize, const unsigned& t_meanDegree, const std::vector<double>& t_initialRates, const double& t_temperature) : m_networkSize(t_networkSize), m_meanDegree(t_meanDegree), m_rates(t_initialRates), m_temperature(t_temperature) {
    //* Generate seed Engine with random initial seed
    m_randomEngineSeedDistribution.param(std::uniform_int_distribution<int>::param_type(0, 10000));
    m_probDistribution.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    m_seedEngine.seed((std::random_device())());

    //* Real real data
    CSV::read("realConfirmed_208.txt", m_realConfirmed);
    m_maxDate = m_realConfirmed.size();

    //* Read initial energy
    const std::string target = "ER,N" + to_stringWithExponent((double)t_networkSize, 1) + ",M" + std::to_string(t_meanDegree) + "/" + "SE" + to_stringWithPrecision(t_initialRates[0], 2) + ",EAI" + to_stringWithPrecision(t_initialRates[1], 2) + ",pA" + to_stringWithPrecision(t_initialRates[2], 2) + ",IQI" + to_stringWithPrecision(t_initialRates[3], 2) + ",AR" + to_stringWithPrecision(t_initialRates[4], 2) + ",QICR" + to_stringWithPrecision(t_initialRates[5], 2) + ",XQX" + to_stringWithPrecision(t_initialRates[6], 2) + ",T" + to_stringWithPrecision(t_initialRates[7], 2);
    std::ifstream energyFile("energyList.txt");
    std::string line;
    while (getline(energyFile, line)){
        if (line.find(target) != line.npos){
            const double energy = std::stod(line.substr(line.find(": ") + 2));
            if (energy < m_energy){
                m_energy = energy;
            }
        }
    }
}

const double SA_KM::m_getEnergy(const std::vector<std::vector<unsigned>>& t_simulatedData) const {
    double RMSE = 0.0;

    //* Get confirmed case from simulated data
    std::vector<unsigned> simulatedConfirmed;
    simulatedConfirmed.reserve(t_simulatedData.size());
    for (const std::vector<unsigned>& simulatedData : t_simulatedData) {
        simulatedConfirmed.emplace_back(simulatedData[7] + simulatedData[8] + simulatedData[10]);
    }

    //* Get RMSE from difference between simulation and real data
    for (unsigned date = 0; date < m_maxDate; ++date) {
        RMSE += std::pow((double)m_realConfirmed[date] - (double)simulatedConfirmed[date], 2.0);
    }
    return std::sqrt(RMSE / (double)m_maxDate);
}

const std::vector<double> SA_KM::m_perturbate(const double& t_perturbationSize = 0.01) {
    std::uniform_int_distribution<int> rateDistribution(0, m_rates.size());

    std::vector<double> newRates = m_rates;
    const int chosenRate = rateDistribution(m_seedEngine);
    m_probDistribution(m_seedEngine) > 0.5 ? newRates[chosenRate] += t_perturbationSize : newRates[chosenRate] -= t_perturbationSize;

    return newRates;
}

void SA_KM::sync_run(const unsigned& t_maxEnsemble, const double& t_deltaT=1e-1) {
    for (unsigned ensemble=0; ensemble < t_maxEnsemble; ++ensemble){
        //* Perturbate rate of KM simulation
        const std::vector<double> newRates = m_perturbate();

        //* Generate random Engine from seed engine
        const int randomEngineSeed = m_randomEngineSeedDistribution(m_seedEngine);
        pcg32 randomEngine(randomEngineSeed);

        //* Generate network
        const Network network = ER::generate(m_networkSize, m_networkSize * m_meanDegree / 2, randomEngine);

        //* Generate KM model and path for data
        KM simulation(network, newRates, randomEngine);
        const std::string dataDirectory = "../data/KM/";
        const std::string networkDirectory = simulation.directoryName("ER", m_networkSize, m_meanDegree);
        CSV::generateDirectory(dataDirectory + networkDirectory);
        const std::string fileName = simulation.fileName(randomEngineSeed);

        //* Run simulation and save data
        double newEnergy = 0.0;
        const std::vector<std::vector<unsigned>> simulatedData = simulation.sync_run(t_deltaT, m_maxDate);
        if (simulatedData.size() < m_maxDate) {
            std::cout << "At random seed " << randomEngineSeed << ", simulated finished before reacing current time.\n";
            continue;
        } else {
            newEnergy = m_getEnergy(simulatedData);
            std::ofstream energyFile("energyList.txt", std::ios_base::app);
            energyFile << networkDirectory << fileName << ": " << newEnergy << "\n";
        }
        CSV::write(dataDirectory + networkDirectory + fileName, simulatedData);

        //* Update rate according to energy
        const double deltaE = newEnergy - m_energy;
        if (m_probDistribution(m_seedEngine) < std::exp(-1.0 * deltaE / m_temperature)){
            m_energy = newEnergy;
            m_rates = newRates;
            std::cout << "Changed rates\n";
        }

    }
}