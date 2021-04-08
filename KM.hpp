#pragma once

#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "../library/CSV.hpp"
#include "../library/Epidemics.hpp"
#include "../library/Networks.hpp"
#include "../library/linearAlgebra.hpp"
#include "../library/pcg_random.hpp"
#include "../library/stringFormat.hpp"

/*
    K-Model simulation
    S + A(I) -> E + A(I) with rate S_E
    E -> A with rate m_EA
    E -> I with rate m_EI
    I -> QI with rate m_IQI
    A(QA) -> R(CR) with rate m_AR
    QI -> CR with rate m_QICR
    X -> QX with rate m_XQX where X=S,E,A,I,R
    QX -> X with constant time tau
*/

const std::map<std::string, int> stateToInt = {{"S", 0}, {"E", 1}, {"A", 2}, {"I", 3}, {"R", 4}, {"QS", 5}, {"QE", 6}, {"QA", 7}, {"QI", 8}, {"QR", 9}, {"CR", 10}};

struct KNode : public Node_Epidemic<unsigned> {
    //* Member variables
    double quarantineTime{0.0};      // Time spend after became quarantined
    unsigned quarantineNeighbor{0};  // Number of QA, QI neighbors
    unsigned infectiousNeighbor{0};

    //* Generator
    KNode() {}
    KNode(const unsigned& t_index) : Node_Epidemic<unsigned>(t_index) {}
    KNode(const unsigned& t_index, const std::string& t_state) : Node_Epidemic<unsigned>(t_index, t_state) {}
};

struct KM {
    //* Member variables
   protected:
    std::vector<KNode> m_nodes;
    std::set<unsigned> m_reactingIndex, m_quarantineIndex;
    double m_currentTime;
    unsigned m_nextDate;
    std::vector<std::vector<unsigned>> m_data;
    double m_SE, m_EA, m_EI, m_IQI, m_AR, m_QICR, m_XQX, m_tau;
    std::vector<unsigned> m_numberOfStates;
    pcg32 m_randomEngine;
    std::uniform_real_distribution<double> m_probabilityDistribution;

    //* Member functions
   protected:
    const unsigned m_getQuarantineNeighbor(const unsigned&) const;
    const unsigned m_getInfectiousNeighbor(const unsigned&) const;
    void m_updateTransitionRate(const unsigned&);
    void m_updateTransitionRates();
    void m_release(const double&);
    void m_syncUpdate(const double&);

   public:
    //* Generator
    KM() {}
    KM(const Network<unsigned>&, const std::vector<double>&, const pcg32&, const int&);
    const bool sync_run(const double&, const unsigned&);
    const double getEnergy(const std::vector<unsigned>&) const;
    void save(const std::string& t_file) const { CSV::write(t_file, m_data); }
};

KM::KM(const Network<unsigned>& t_network, const std::vector<double>& t_rates, const pcg32& t_randomEngine, const int& t_seedSize = 1) : m_randomEngine(t_randomEngine) {
    m_probabilityDistribution.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));

    //* Set Network
    const unsigned networkSize = t_network.size;
    m_nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        KNode node(index, "S");
        node.neighbors = t_network.adjacency[index];
        m_nodes.emplace_back(node);
    }

    //* Set rates
    m_SE = t_rates[0];
    m_EA = t_rates[1] * t_rates[2];
    m_EI = t_rates[1] * (1 - t_rates[2]);
    m_IQI = t_rates[3];
    m_AR = t_rates[4];
    m_QICR = t_rates[5];
    m_XQX = t_rates[6];
    m_tau = t_rates[7];

    //* Set initial seed and number of states
    m_numberOfStates.assign(stateToInt.size(), 0);
    m_numberOfStates[0] = networkSize - t_seedSize;
    m_numberOfStates[1] = t_seedSize;
    for (int seed = 0; seed < t_seedSize; ++seed) {
        m_nodes[seed].state = "E";
        m_reactingIndex.emplace(seed);
    }
    m_data.reserve(500);
    m_data.emplace_back(m_numberOfStates);

    //* Set transition rates of nodes and time
    m_currentTime = 0.0;
    m_nextDate = 1;
    m_updateTransitionRates();
}

const unsigned KM::m_getQuarantineNeighbor(const unsigned& t_index) const {
    unsigned quarantineNeighbor = 0;
    for (const unsigned& neighbor : m_nodes[t_index].neighbors) {
        if (m_nodes[neighbor].state == "QA" || m_nodes[neighbor].state == "QI") {
            ++quarantineNeighbor;
        }
    }
    return quarantineNeighbor;
}

const unsigned KM::m_getInfectiousNeighbor(const unsigned& t_index) const {
    unsigned infectiousNeighbor = 0;
    for (const unsigned& neighbor : m_nodes[t_index].neighbors) {
        if (m_nodes[neighbor].state == "A" || m_nodes[neighbor].state == "I") {
            ++infectiousNeighbor;
        }
    }
    return infectiousNeighbor;
}

void KM::m_updateTransitionRate(const unsigned& t_index) {
    //* Update number of neighbor QA, QI
    const unsigned quarantineNeighbor = m_getQuarantineNeighbor(t_index);
    m_nodes[t_index].quarantineNeighbor = quarantineNeighbor;
    const int intState = stateToInt.at(m_nodes[t_index].state);
    switch (intState) {
        //* S process
        case 0: {
            const unsigned infectiousNeighbor = m_getInfectiousNeighbor(t_index);
            m_nodes[t_index].infectiousNeighbor = infectiousNeighbor;
            m_nodes[t_index].transitionRate = m_SE * infectiousNeighbor + m_XQX * quarantineNeighbor;
            break;
        }
        //* E process
        case 1: {
            m_nodes[t_index].transitionRate = m_EA + m_EI + m_XQX * quarantineNeighbor;
            break;
        }
        //* A process
        case 2: {
            m_nodes[t_index].transitionRate = m_AR + m_XQX * quarantineNeighbor;
            break;
        }
        //* I process
        case 3: {
            m_nodes[t_index].transitionRate = m_IQI + m_XQX * quarantineNeighbor;
            break;
        }
        //* R process
        case 4: {
            m_nodes[t_index].transitionRate = m_XQX * quarantineNeighbor;
            break;
        }
        //* QE process
        case 6: {
            m_nodes[t_index].transitionRate = m_EA + m_EI;
            break;
        }
        //* QA process
        case 7: {
            m_nodes[t_index].transitionRate = m_AR;
            break;
        }
        //* QI process
        case 8: {
            m_nodes[t_index].transitionRate = m_QICR;
            break;
        }
        //* QS, QR, CR process
        default: {
            m_nodes[t_index].transitionRate = 0.0;
            break;
        }
    }
}

void KM::m_updateTransitionRates() {
    for (const unsigned& index : m_reactingIndex) {
        m_updateTransitionRate(index);
    }
}

void KM::m_release(const double& t_deltaT) {
    const std::set<unsigned> quarantineIndex = m_quarantineIndex;
    for (const unsigned& index : quarantineIndex) {
        m_nodes[index].quarantineTime += t_deltaT;
        const int intState = stateToInt.at(m_nodes[index].state);
        switch (intState) {
            //* QS process
            case 5: {
                if (m_nodes[index].quarantineTime > m_tau) {
                    m_nodes[index].state = "S";
                    --m_numberOfStates[5];
                    ++m_numberOfStates[0];
                    m_nodes[index].quarantineTime = 0.0;
                    m_quarantineIndex.erase(index);
                    m_updateTransitionRate(index);
                    if (m_nodes[index].transitionRate) {
                        m_reactingIndex.emplace(index);
                    }
                }
                break;
            }
            //* QE process
            case 6: {
                if (m_nodes[index].quarantineTime > m_tau) {
                    m_nodes[index].state = "E";
                    --m_numberOfStates[6];
                    ++m_numberOfStates[1];
                    m_nodes[index].quarantineTime = 0.0;
                    m_quarantineIndex.erase(index);
                    m_updateTransitionRate(index);
                }
                break;
            }
            //* QR process
            case 9: {
                if (m_nodes[index].quarantineTime > m_tau) {
                    m_nodes[index].state = "R";
                    --m_numberOfStates[9];
                    ++m_numberOfStates[4];
                    m_nodes[index].quarantineTime = 0.0;
                    m_quarantineIndex.erase(index);
                    m_updateTransitionRate(index);
                    if (m_nodes[index].transitionRate) {
                        m_reactingIndex.emplace(index);
                    }
                }
                break;
            }
            //* QA, QI process
            default: {
                break;
            }
        }
    }
}

void KM::m_syncUpdate(const double& t_deltaT) {
    //* Update time and check to store data
    bool storeData = false;
    m_currentTime += t_deltaT;
    if (m_currentTime >= m_nextDate){
        ++m_nextDate;
        storeData = true;
    }

    //* Release from quarantine
    m_release(t_deltaT);

    //* Do reactions according to each transition rate and add E,A,I,QE,QA,QI into reacting nodes
    std::set<unsigned> newReactingIndex = m_reactingIndex;
    for (const unsigned& index : m_reactingIndex) {
        const int intState = stateToInt.at(m_nodes[index].state);
        const double transitionRate = m_nodes[index].transitionRate;
        const double transitionProb = 1.0 - std::exp(-1.0 * transitionRate * t_deltaT);
        const double prob = m_probabilityDistribution(m_randomEngine);

        switch (intState) {
            //* S process
            case 0: {
                if (prob <= transitionProb) {
                    const double p = m_probabilityDistribution(m_randomEngine);
                    //* S -> QS
                    if (p * transitionRate < m_nodes[index].quarantineNeighbor * m_XQX) {
                        m_nodes[index].state = "QS";
                        --m_numberOfStates[0];
                        ++m_numberOfStates[5];
                        m_nodes[index].quarantineTime = 0.0;
                        m_quarantineIndex.emplace(index);
                        newReactingIndex.erase(index);
                    }
                    //* S -> E
                    else {
                        m_nodes[index].state = "E";
                        --m_numberOfStates[0];
                        ++m_numberOfStates[1];
                    }
                }
                //* S -> S
                else {
                    newReactingIndex.erase(index);
                }
                break;
            }
            //* E process
            case 1: {
                if (prob <= transitionProb) {
                    const double p = m_probabilityDistribution(m_randomEngine);
                    //* E -> A
                    if (p * transitionRate < m_EA) {
                        m_nodes[index].state = "A";
                        --m_numberOfStates[1];
                        ++m_numberOfStates[2];
                    }
                    //* E -> I
                    else if (p * transitionRate < m_EA + m_EI) {
                        m_nodes[index].state = "I";
                        --m_numberOfStates[1];
                        ++m_numberOfStates[3];
                    }
                    //* E -> QE
                    else {
                        m_nodes[index].state = "QE";
                        --m_numberOfStates[1];
                        ++m_numberOfStates[6];
                        m_nodes[index].quarantineTime = 0.0;
                        m_quarantineIndex.emplace(index);
                    }
                }
                //* E -> E
                break;
            }
            //* A process
            case 2: {
                if (prob <= transitionProb) {
                    const double p = m_probabilityDistribution(m_randomEngine);
                    //* A -> R
                    if (p * transitionRate < m_AR) {
                        m_nodes[index].state = "R";
                        --m_numberOfStates[2];
                        ++m_numberOfStates[4];
                        newReactingIndex.erase(index);
                    }
                    //* A->QA
                    else {
                        m_nodes[index].state = "QA";
                        --m_numberOfStates[2];
                        ++m_numberOfStates[7];
                    }
                }
                //* A -> A
                break;
            }
            //* I process
            case 3: {
                //* I -> QI
                if (prob <= transitionProb) {
                    m_nodes[index].state = "QI";
                    --m_numberOfStates[3];
                    ++m_numberOfStates[8];
                }
                //* I -> I
                break;
            }
            //* R process
            case 4: {
                //* R -> QR
                if (prob <= transitionProb) {
                    m_nodes[index].state = "QR";
                    --m_numberOfStates[4];
                    ++m_numberOfStates[9];
                    m_nodes[index].quarantineTime = 0.0;
                    m_quarantineIndex.emplace(index);
                }
                //* R -> R
                newReactingIndex.erase(index);
                break;
            }
            //* QE process
            case 6: {
                if (prob <= transitionProb) {
                    const double p = m_probabilityDistribution(m_randomEngine);
                    //* QE -> QA
                    if (p * transitionRate < m_EA) {
                        m_nodes[index].state = "QA";
                        --m_numberOfStates[6];
                        ++m_numberOfStates[7];
                    }
                    //* QE -> QI
                    else {
                        m_nodes[index].state = "QI";
                        --m_numberOfStates[6];
                        ++m_numberOfStates[8];
                    }
                }
                //* QE -> QE
                break;
            }
            //* QA process
            case 7: {
                //* QA -> CR
                if (prob <= transitionProb) {
                    m_nodes[index].state = "CR";
                    --m_numberOfStates[7];
                    ++m_numberOfStates[10];
                    newReactingIndex.erase(index);
                }
                //* QA -> QA
                break;
            }
            //* QI process
            case 8: {
                //* QI -> CR
                if (prob <= transitionProb) {
                    m_nodes[index].state = "CR";
                    --m_numberOfStates[8];
                    ++m_numberOfStates[10];
                    newReactingIndex.erase(index);
                }
                //* QI -> QI
                break;
            }
            //* QS, QR, CR process
            default: {
                break;
            }
        }
    }

    m_reactingIndex = newReactingIndex;
    //* Add neighbors to reacting index
    for (const unsigned& index : newReactingIndex) {
        //* Add neighbor S of A,I node into reacting nodes (spreading)
        if (m_nodes[index].state == "A" || m_nodes[index].state == "I") {
            for (const unsigned& neighbor : m_nodes[index].neighbors) {
                if (m_nodes[neighbor].state == "S") {
                    m_reactingIndex.emplace(neighbor);
                }
            }
        }
        //* Add neighbor S,R of QA,QI node into reacting nodes (quarantine)
        else if (m_nodes[index].state == "QA" || m_nodes[index].state == "QI") {
            for (const unsigned& neighbor : m_nodes[index].neighbors) {
                if (m_nodes[neighbor].state == "S" || m_nodes[neighbor].state == "R") {
                    m_reactingIndex.emplace(neighbor);
                }
            }
        }
    }

    //* Update time and Transition rate
    if (storeData){
        m_data.emplace_back(m_numberOfStates);
    }

    //* Update transition rate of every reacting nodes
    m_updateTransitionRates();
}

const bool KM::sync_run(const double& t_deltaT, const unsigned& t_maxDate) {
    while (m_nextDate < t_maxDate) {
        m_syncUpdate(t_deltaT);
        if (m_reactingIndex.empty()) {
            m_data.emplace_back(m_numberOfStates);
            return false;
        }
    }
    return true;
}

const double KM::getEnergy(const std::vector<unsigned>& t_realConfirmed) const {
    const unsigned maxDate = t_realConfirmed.size();
    if (m_data.size() != maxDate) {
        std::cout << "In get energy function, length of data does not corresponds to max date\n";
        exit(1);
    }

    //* Get RMSE from difference between simulation and real data
    double RMSE = 0.0;
    for (unsigned date = 0; date < maxDate; ++date) {
        RMSE += std::pow((double)t_realConfirmed[date] - (double)m_data[date][7] - (double)m_data[date][8] - (double)m_data[date][10], 2.0);
    }
    RMSE = std::sqrt(RMSE / (double)maxDate);
    return RMSE;
}
