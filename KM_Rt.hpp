#pragma once

#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "../library/Epidemics.hpp"
#include "../library/Networks.hpp"
#include "../library/pcg_random.hpp"
#include "../library/stringFormat.hpp"

/*
    R(t) : Number of second infection from node who was infected (became E) at time [t, t+1)
*/

const std::map<std::string, int> stateToInt = {{"S", 0}, {"E", 1}, {"A", 2}, {"I", 3}, {"R", 4}, {"QS", 5}, {"QE", 6}, {"QA", 7}, {"QI", 8}, {"QR", 9}, {"CR", 10}};

struct KNode_Rt : public Node_Epidemic {
    //* Member variables
    double quarantineTime{0.0};      // Time spend after became quarantined
    unsigned quarantineNeighbor{0};  // Number of QA, QI neighbors
    unsigned infectiousNeighbor{0};
    unsigned infectedDate{0};
    double avergeInfect{0.0};
    std::set<unsigned> infectiousNeighborIndex;

    KNode_Rt() {}
    KNode_Rt(const unsigned& t_index) : Node_Epidemic(t_index) {}
    KNode_Rt(const unsigned& t_index, const std::string& t_state) : Node_Epidemic(t_index, t_state) {}
};

struct KM_Rt {
    //* Member variables
   protected:
    std::vector<KNode_Rt> m_nodes;
    std::vector<std::set<unsigned>> m_currentAdjacency;
    std::set<unsigned> m_reactingIndex, m_quarantineIndex;
    unsigned m_seedSize;
    double m_currentTime;
    unsigned m_nextDate;
    std::vector<std::vector<unsigned>> m_data;
    double m_SE, m_EA, m_EI, m_IQI, m_AR, m_QICR, m_XQX, m_tau;
    std::vector<unsigned> m_numberOfStates;                  //* S, E, A, I, R, QS, QE, QA, QI, QR, CR
    std::vector<std::pair<unsigned, double>> m_avgLinkSize;  //* Number of intervals, sum of every link size
    pcg32 m_randomEngine;
    std::uniform_real_distribution<double> m_probabilityDistribution;

    //* member functions
   protected:
    const unsigned m_getQuarantineNeighbor(const unsigned&) const;
    const unsigned m_getInfectiousNeighbor(const unsigned&);
    void m_updateTransitionRate(const unsigned&);
    void m_updateTransitionRates();
    void m_isolateNode(const unsigned&);
    void m_restoreNode(const unsigned&);
    void m_release(const double&);
    void m_syncUpdate(const double&);
    const unsigned m_getLinkSize() const;

   public:
    //* Generator
    KM_Rt() {}
    KM_Rt(const Network&, const std::vector<double>&, const pcg32&, const int&);
    const bool sync_run(const double&, const unsigned&);
    void save(const std::string&) const;
};

KM_Rt::KM_Rt(const Network& t_network, const std::vector<double>& t_rates, const pcg32& t_randomEngine, const int& t_seedSize = 1) : m_randomEngine(t_randomEngine), m_seedSize(t_seedSize) {
    m_probabilityDistribution.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));

    //* Set Network
    const unsigned networkSize = t_network.m_size;
    m_currentAdjacency = t_network.m_adjacency;
    m_nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        KNode_Rt node(index, "S");
        node.m_neighbors = t_network.m_adjacency[index];
        m_nodes.emplace_back(node);
    }
    m_avgLinkSize.reserve(500);
    m_avgLinkSize.emplace_back(std::make_pair<unsigned, double>(1, (double)t_network.m_linkSize));

    //* Set rates
    m_SE = t_rates[0];
    m_EA = t_rates[1] * t_rates[2];
    m_EI = t_rates[1] * (1 - t_rates[2]);
    m_IQI = t_rates[3];
    m_AR = t_rates[4];
    m_QICR = t_rates[5];
    m_XQX = t_rates[6];
    m_tau = t_rates[7];

    //* Set initial seed
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

const unsigned KM_Rt::m_getLinkSize() const {
    unsigned linkSize = 0;
    for (const std::set<unsigned>& adjacent : m_currentAdjacency) {
        linkSize += adjacent.size();
    }
    return linkSize / 2;
}

const unsigned KM_Rt::m_getQuarantineNeighbor(const unsigned& t_index) const {
    unsigned quarantineNeighbor = 0;
    for (const unsigned& neighbor : m_nodes[t_index].m_neighbors) {
        if (m_nodes[neighbor].state == "QA" || m_nodes[neighbor].state == "QI") {
            ++quarantineNeighbor;
        }
    }
    return quarantineNeighbor;
}

const unsigned KM_Rt::m_getInfectiousNeighbor(const unsigned& t_index) {
    unsigned infectiousNeighbor = 0;
    m_nodes[t_index].infectiousNeighborIndex.clear();
    for (const unsigned& neighbor : m_nodes[t_index].m_neighbors) {
        if (m_nodes[neighbor].state == "A" || m_nodes[neighbor].state == "I") {
            ++infectiousNeighbor;
            m_nodes[t_index].infectiousNeighborIndex.emplace(neighbor);
        }
    }
    return infectiousNeighbor;
}

void KM_Rt::m_updateTransitionRate(const unsigned& t_index) {
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

void KM_Rt::m_updateTransitionRates() {
    for (const unsigned& index : m_reactingIndex) {
        m_updateTransitionRate(index);
    }
}

void KM_Rt::m_isolateNode(const unsigned& t_index) {
    for (const unsigned& neighbor : m_nodes[t_index].m_neighbors) {
        m_currentAdjacency[t_index].erase(neighbor);
        m_currentAdjacency[neighbor].erase(t_index);
    }
}

void KM_Rt::m_restoreNode(const unsigned& t_index) {
    for (const unsigned& neighbor : m_nodes[t_index].m_neighbors) {
        if (m_nodes[neighbor].state.find("Q") == m_nodes[neighbor].state.npos) {
            m_currentAdjacency[t_index].emplace(neighbor);
            m_currentAdjacency[neighbor].emplace(t_index);
        }
    }
}

void KM_Rt::m_release(const double& t_deltaT) {
    const std::set<unsigned> quarantineIndex = m_quarantineIndex;
    for (const unsigned& index : quarantineIndex) {
        m_nodes[index].quarantineTime += t_deltaT;
        const int intState = stateToInt.at(m_nodes[index].state);
        switch (intState) {
            //* QS -> S
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
                    m_restoreNode(index);
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
                    m_restoreNode(index);
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
                    m_restoreNode(index);
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

void KM_Rt::m_syncUpdate(const double& t_deltaT) {
    //* Update time and check to store data
    bool storeData = false;
    m_currentTime += t_deltaT;
    if (m_currentTime >= m_nextDate) {
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
                        m_isolateNode(index);
                    }
                    //* S -> E
                    else {
                        m_nodes[index].state = "E";
                        --m_numberOfStates[0];
                        ++m_numberOfStates[1];
                        //? Changed for reproduction number
                        m_nodes[index].infectedDate = m_nextDate - 1;
                        for (const unsigned& neighbor : m_nodes[index].infectiousNeighborIndex) {
                            m_nodes[neighbor].avergeInfect += 1.0 / (double)m_nodes[index].infectiousNeighbor;
                        }
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
                        m_isolateNode(index);
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
                        m_isolateNode(index);
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
                    m_isolateNode(index);
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
                    m_isolateNode(index);
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
            for (const unsigned& neighbor : m_nodes[index].m_neighbors) {
                if (m_nodes[neighbor].state == "S") {
                    m_reactingIndex.emplace(neighbor);
                }
            }
        }
        //* Add neighbor S,R of QA,QI node into reacting nodes (quarantine)
        else if (m_nodes[index].state == "QA" || m_nodes[index].state == "QI") {
            for (const unsigned& neighbor : m_nodes[index].m_neighbors) {
                if (m_nodes[neighbor].state == "S" || m_nodes[neighbor].state == "R") {
                    m_reactingIndex.emplace(neighbor);
                }
            }
        }
    }

    //* Store data if neccessary and update average link size
    //? Changed for reproduction number
    if (storeData) {
        m_data.emplace_back(m_numberOfStates);
        m_avgLinkSize.emplace_back(std::make_pair<unsigned, double>(1, (double)m_getLinkSize()));
    } else {
        ++m_avgLinkSize[m_nextDate - 1].first;
        m_avgLinkSize[m_nextDate - 1].second += (double)m_getLinkSize();
    }

    //* Update transition rate of every reacting nodes
    m_updateTransitionRates();
}

const bool KM_Rt::sync_run(const double& t_deltaT, const unsigned& t_maxDate) {
    while (m_nextDate < t_maxDate) {
        m_syncUpdate(t_deltaT);
        if (m_reactingIndex.empty()) {
            m_data.emplace_back(m_numberOfStates);
            return false;
        }
    }
    return true;
}

void KM_Rt::save(const std::string& t_file) const {
    std::vector<std::vector<double>> totalData;
    totalData.reserve(m_data.size());

    //* Follow nodes that have not counted as reproduction number yet
    std::set<unsigned> remainedNodesIndex;
    for (unsigned index = 0; index < m_nodes.size(); ++index) {
        remainedNodesIndex.emplace_hint(remainedNodesIndex.end(), index);
    }

    //* From day1, consider every nodes
    for (unsigned date = 0; date < m_data.size(); ++date) {
        //* Save numer of states
        std::vector<double> temp(m_data[date].begin(), m_data[date].end());

        //* Save average link size
        const double avgLink = m_avgLinkSize[date].second / (double)m_avgLinkSize[date].first;
        temp.emplace_back(avgLink);

        //* Save reproduction number
        double reproduction = 0.0;
        unsigned num = 0;
        const std::set<unsigned> temp_remainedNodesIndex = remainedNodesIndex;
        for (const unsigned& index : temp_remainedNodesIndex) {
            if (m_nodes[index].infectedDate == date) {
                ++num;
                reproduction += m_nodes[index].avergeInfect;
                remainedNodesIndex.erase(index);
            }
        }
        if (date == 0){
            num = m_seedSize;
        }
        reproduction = num > 0 ? reproduction / num : 0.0;
        temp.emplace_back(reproduction);

        //* Save to total data
        totalData.emplace_back(temp);
    }

    CSV::write(t_file, totalData);
}
