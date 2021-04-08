#pragma once

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "../library/CSV.hpp"
#include "../library/Networks.hpp"
#include "../library/pcg_random.hpp"

#include "KM.hpp"
#include "fileName.hpp"

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
    const std::vector<double> m_perturbate(const double&);

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
    CSV::read("../data/COVID/realData/confirmed_209.txt", m_realConfirmed);
    m_maxDate = m_realConfirmed.size();

    //* Read initial energy
    const std::string target = "ER,N" + to_stringWithExponent((double)t_networkSize, 1) + ",M" + std::to_string(t_meanDegree) + "/" + "SE" + to_stringWithPrecision(t_initialRates[0], 2) + ",EAI" + to_stringWithPrecision(t_initialRates[1], 2) + ",pA" + to_stringWithPrecision(t_initialRates[2], 2) + ",IQI" + to_stringWithPrecision(t_initialRates[3], 2) + ",AR" + to_stringWithPrecision(t_initialRates[4], 2) + ",QICR" + to_stringWithPrecision(t_initialRates[5], 2) + ",XQX" + to_stringWithPrecision(t_initialRates[6], 2) + ",T" + to_stringWithPrecision(t_initialRates[7], 2);
    std::ifstream energyFile("energyList.txt");
    std::string line;
    while (getline(energyFile, line)) {
        if (line.find(target) != line.npos) {
            const double energy = std::stod(line.substr(line.find(": ") + 2));
            if (energy < m_energy) {
                m_energy = energy;
            }
        }
    }
}

const std::vector<double> SA_KM::m_perturbate(const double& t_perturbationSize = 0.001) {
    std::vector<double> newRates = m_rates;

    //* Perturbate only rate[0]: SE or rate[6]: XQX
    const int chosenRate = m_probDistribution(m_seedEngine) < 0.5 ? 0 : 6;
    m_probDistribution(m_seedEngine) > 0.5 ? newRates[chosenRate] += t_perturbationSize : newRates[chosenRate] -= t_perturbationSize;

    return newRates;
}

void SA_KM::sync_run(const unsigned& t_maxEnsemble, const double& t_deltaT = 1e-1) {
    for (unsigned temperatureIndex=0; temperatureIndex<10; ++temperatureIndex){
        for (unsigned ensemble = 0; ensemble < t_maxEnsemble; ++ensemble) {
            //* Perturbate rate of KM simulation
            const std::vector<double> newRates = m_perturbate(0.001);

            //* Generate random Engine from seed engine
            const int randomEngineSeed = m_randomEngineSeedDistribution(m_seedEngine);
            pcg32 randomEngine(randomEngineSeed);

            //* Generate network
            const Network<unsigned> network = ER::generate(m_networkSize, (unsigned long long)m_networkSize * m_meanDegree / 2, randomEngine);

            //* Generate KM model and path for data
            KM simulation(network, newRates, randomEngine);
            const std::string dataDirectory = "../data/KM/";
            const std::string networkDirectory = networkName("ER", m_networkSize, m_meanDegree);
            const std::string fileName = rateName(newRates, randomEngineSeed);
            CSV::generateDirectory(dataDirectory + networkDirectory);

            //* Run simulation and save data
            double newEnergy = 0.0;
            const bool finished = simulation.sync_run(t_deltaT, m_maxDate);
            if (finished) {
                newEnergy = simulation.getEnergy(m_realConfirmed);
            } else {
                // std::cout << "At random seed " << randomEngineSeed << ", simulated finished before reacing current time.\n";
                continue;
            }

            //* Update rate according to energy
            const double deltaE = newEnergy - m_energy;
            if (m_probDistribution(m_seedEngine) < std::exp(-1.0 * deltaE / m_temperature)) {
                m_energy = newEnergy;
                m_rates = newRates;
                std::ofstream energyFile("energyList.txt", std::ios_base::app);
                energyFile << networkDirectory << fileName << ": " << newEnergy << "\n";
                simulation.save(dataDirectory + networkDirectory + fileName);
                // std::cout << "Changed rates\n";
            }
        }
        m_temperature /= 2;
    }
}