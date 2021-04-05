#pragma once

#include <iostream>
#include <map>
#include <random>  //uniform distribution
#include <set>
#include <string>
#include <vector>

#include "../library/CSV.hpp"
#include "../library/Epidemics.hpp"
#include "../library/Networks.hpp"
#include "../library/pcg_random.hpp"
#include "../library/stringFormat.hpp"

/*
    K-Model simulation
    S + A(I) -> E + A(I) with rate S_E
    E -> A with rate E_A
    E -> I with rate E_I
    I -> QI with rate I_QI
    A(QA) -> R(CR) with rate A_R
    QI -> CR with rate QI_CR
    X -> QX with rate X_QX where X=S,E,A,I,R
    QX -> X with constant time tau
*/

struct KNode : public Node_Epidemic {
    //* Member variables
    double m_quarantineTime{0.0};      // Time spend after became quarantined
    unsigned m_quarantineNeighbor{0};  // Number of QA, QI neighbors

    //* Generator
    KNode() {}
    KNode(const unsigned& t_index) : Node_Epidemic(t_index) {}
    KNode(const unsigned& t_index, const std::string& t_state) : Node_Epidemic(t_index, t_state) {}
};

namespace KM {
//* Save variables
bool deletion;
std::string directory;

//* Rate parameters
double S_E,E_AI, pA, E_A, E_I, I_QI, A_R, QI_CR, X_QX, tau;

//* Random variables
int randomEngineSeed;
pcg32 randomEngine;
std::uniform_real_distribution<double> probabilityDistribution(0, 1);

//* Parameters used for GA model
unsigned networkSize;
const std::map<std::string, int> stateToInt = {{"S", 0}, {"E", 1}, {"A", 2}, {"I", 3}, {"R", 4}, {"QS", 5}, {"QE", 6}, {"QA", 7}, {"QI", 8}, {"QR", 9}, {"CR", 10}};
unsigned seedSize;
std::vector<KNode> nodes;
std::set<unsigned> reactingIndex;
std::set<unsigned> quarantineIndex;  //* Track quarantine time of QS, QE, QR for release
double currentTime;
unsigned numS, numE, numA, numI, numR, numQS, numQE, numQA, numQI, numQR, numCR;

//* File name convention
const std::string fileName() {
    const std::string fileName = "SE" + to_stringWithPrecision(S_E, 2) + ",EAI" + to_stringWithPrecision(E_AI, 2) + ",pA" + to_stringWithPrecision(pA, 2) + ",IQI" + to_stringWithPrecision(I_QI, 2) + ",AR" + to_stringWithPrecision(A_R, 2) + ",QICR" + to_stringWithPrecision(QI_CR, 2) + ",XQX" + to_stringWithPrecision(X_QX, 2) + ",T" + to_stringWithPrecision(tau, 2);

    return randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(randomEngineSeed) + ".txt";
}

//* Update transition rate of single node
void updateTransitionRate(const unsigned& t_index) {
    //* Update number of neighbor QA, QI
    unsigned quarantineNeighbor = 0;
    for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
        if (nodes[neighbor].m_state == "QA" || nodes[neighbor].m_state == "QI") {
            ++quarantineNeighbor;
        }
    }
    nodes[t_index].m_quarantineNeighbor = quarantineNeighbor;
    const int intState = stateToInt.at(nodes[t_index].m_state);
    switch (intState) {
        //* S process
        case 0: {
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
                if (nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].m_transitionRate = S_E * infectiousNeighbor + X_QX * quarantineNeighbor;
            break;
        }
        //* E process
        case 1: {
            nodes[t_index].m_transitionRate = E_A + E_I + X_QX * quarantineNeighbor;
            break;
        }
        //* A process
        case 2: {
            nodes[t_index].m_transitionRate = A_R + X_QX * quarantineNeighbor;
            break;
        }
        //* I process
        case 3: {
            nodes[t_index].m_transitionRate = I_QI + X_QX * quarantineNeighbor;
            break;
        }
        //* R process
        case 4: {
            nodes[t_index].m_transitionRate = X_QX * quarantineNeighbor;
            break;
        }
        //* QE process
        case 6: {
            nodes[t_index].m_transitionRate = E_A + E_I;
            break;
        }
        //* QA process
        case 7: {
            nodes[t_index].m_transitionRate = A_R;
            break;
        }
        //* QI process
        case 8: {
            nodes[t_index].m_transitionRate = QI_CR;
            break;
        }
        //* QS, QR, CR process
        default: {
            nodes[t_index].m_transitionRate = 0.0;
            break;
        }
    }
}

//* Update transition rate of all reacting nodes
void updateTransitionRate() {
    for (const unsigned& index : reactingIndex) {
        updateTransitionRate(index);
    }
}

//* Initialize model
void initialize(const Network& t_network, const pcg32& t_randomEngine) {
    //! Set nodes
    networkSize = t_network.m_size;
    nodes.clear();
    nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        KNode node(index, "S");
        node.m_neighbors = t_network.m_adjacency[index];
        nodes.emplace_back(node);
    }

    //! Set random Engine
    randomEngine = t_randomEngine;

    //! Initialize model
    currentTime = 0.0;
    numS = networkSize - seedSize;
    numE = seedSize;
    numA = 0;
    numI = 0;
    numR = 0;
    numQS = 0;
    numQE = 0;
    numQA = 0;
    numQI = 0;
    numQR = 0;
    numCR = 0;
    reactingIndex.clear();
    quarantineIndex.clear();
    for (unsigned index = 0; index < seedSize; ++index) {
        nodes[index].m_state = "E";
        reactingIndex.emplace(index);
    }
    updateTransitionRate();
}  //* End of function KM::initialize

//* Release from quarantine
void release(const double& t_deltaT) {
    const std::set<unsigned> temp_quarantineIndex = quarantineIndex;
    for (const unsigned& index : temp_quarantineIndex) {
        nodes[index].m_quarantineTime += t_deltaT;
        const int intState = stateToInt.at(nodes[index].m_state);
        switch (intState) {
            //* QS process
            case 5: {
                if (nodes[index].m_quarantineTime > tau) {
                    nodes[index].m_state = "S";
                    --numQS;
                    ++numS;
                    nodes[index].m_quarantineTime = 0.0;
                    quarantineIndex.erase(index);
                    updateTransitionRate(index);
                    if (nodes[index].m_transitionRate) {
                        reactingIndex.emplace(index);
                    }
                }
                break;
            }
            //* QE process
            case 6: {
                if (nodes[index].m_quarantineTime > tau) {
                    nodes[index].m_state = "E";
                    --numQE;
                    ++numE;
                    nodes[index].m_quarantineTime = 0.0;
                    quarantineIndex.erase(index);
                    updateTransitionRate(index);
                }
                break;
            }
            //* QR process
            case 9: {
                if (nodes[index].m_quarantineTime > tau) {
                    nodes[index].m_state = "R";
                    --numQR;
                    ++numR;
                    nodes[index].m_quarantineTime = 0.0;
                    quarantineIndex.erase(index);
                    updateTransitionRate(index);
                    if (nodes[index].m_transitionRate) {
                        reactingIndex.emplace(index);
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
}  //* End of function KM::release

//* Update one step for every nodes
void syncUpdate(const double& t_deltaT) {
    //* Update current time
    // const double deltaT = t_deltaT;
    currentTime += t_deltaT;

    //* Release from quarantine
    release(t_deltaT);

    //* Do reactions according to each transition rate and add E,A,I,QE,QA,QI into reacting nodes
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : reactingIndex) {
        const int intState = stateToInt.at(nodes[index].m_state);
        const double transitionRate = nodes[index].m_transitionRate;
        const double transitionProb = 1.0 - std::exp(-1.0 * transitionRate * t_deltaT);
        const double prob = probabilityDistribution(randomEngine);
        switch (intState) {
            //* S process
            case 0: {
                if (prob <= transitionProb) {
                    const double p = probabilityDistribution(randomEngine);
                    //* S -> QS
                    if (p * transitionRate < nodes[index].m_quarantineNeighbor * X_QX) {
                        nodes[index].m_state = "QS";
                        --numS;
                        ++numQS;
                        nodes[index].m_transitionRate = 0.0;
                        nodes[index].m_quarantineTime = 0.0;
                        quarantineIndex.emplace(index);
                    }
                    //* S -> E
                    else {
                        nodes[index].m_state = "E";
                        --numS;
                        ++numE;
                        newReactingIndex.emplace(index);
                    }
                }
                break;
            }
            //* E process
            case 1: {
                if (prob <= transitionProb) {
                    const double p = probabilityDistribution(randomEngine);
                    //* E -> A
                    if (p * transitionRate < E_A) {
                        nodes[index].m_state = "A";
                        --numE;
                        ++numA;
                    }
                    //* E -> I
                    else if (p * transitionRate < E_A + E_I) {
                        nodes[index].m_state = "I";
                        --numE;
                        ++numI;
                    }
                    //* E -> QE
                    else {
                        nodes[index].m_state = "QE";
                        --numE;
                        ++numQE;
                        nodes[index].m_quarantineTime = 0.0;
                        quarantineIndex.emplace(index);
                    }
                }
                newReactingIndex.emplace(index);
                break;
            }
            //* A process
            case 2: {
                if (prob <= transitionProb) {
                    const double p = probabilityDistribution(randomEngine);
                    //* A -> R
                    if (p * transitionRate < A_R) {
                        nodes[index].m_state = "R";
                        --numA;
                        ++numR;
                    }
                    //* A->QA
                    else {
                        nodes[index].m_state = "QA";
                        --numA;
                        ++numQA;
                        newReactingIndex.emplace(index);
                    }
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* I process
            case 3: {
                if (prob <= transitionProb) {
                    nodes[index].m_state = "QI";
                    --numI;
                    ++numQI;
                }
                newReactingIndex.emplace(index);
                break;
            }
            //* R process
            case 4: {
                if (prob <= transitionProb) {
                    nodes[index].m_state = "QR";
                    --numR;
                    ++numQR;
                    nodes[index].m_quarantineTime = 0.0;
                    quarantineIndex.emplace(index);
                }
                break;
            }
            //* QE process
            case 6: {
                if (prob <= transitionProb) {
                    const double p = probabilityDistribution(randomEngine);
                    //* QE -> QA
                    if (p * transitionRate < E_A) {
                        nodes[index].m_state = "QA";
                        --numQE;
                        ++numQA;
                    }
                    //* QE -> QI
                    else {
                        nodes[index].m_state = "QI";
                        --numQE;
                        ++numQI;
                    }
                }
                newReactingIndex.emplace(index);
                break;
            }
            //* QA process
            case 7: {
                if (prob <= transitionProb) {
                    nodes[index].m_state = "CR";
                    --numQA;
                    ++numCR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* QI process
            case 8: {
                if (prob <= transitionProb) {
                    nodes[index].m_state = "CR";
                    --numQI;
                    ++numCR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* QS, QR, CR process
            default: {
                break;
            }
        }
    }

    reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        //* Add neighbor S of A,I node into reacting nodes (spreading)
        if (nodes[index].m_state == "A" || nodes[index].m_state == "I") {
            for (const unsigned& neighbor : nodes[index].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    reactingIndex.emplace(neighbor);
                }
            }
        }
        //* Add neighbor S,R of QA,QI node into reacting nodes (quarantine)
        else if (nodes[index].m_state == "QA" || nodes[index].m_state == "QI") {
            for (const unsigned& neighbor : nodes[index].m_neighbors) {
                if (nodes[neighbor].m_state == "S" || nodes[neighbor].m_state == "R") {
                    reactingIndex.emplace(neighbor);
                }
            }
        }
    }

    //* Update time and Transition rate
    updateTransitionRate();

}  //* End of function KM::syncUpdate

double syncRun(const double& t_deltaT) {
    const std::string writeFileName = directory + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS << "," << numE << "," << numA << "," << numI << "," << numR << "," << numQS << "," << numQE << "," << numQA << "," << numQI << "," << numQR << "," << numCR << "\n";
    while (!reactingIndex.empty()) {
        syncUpdate(t_deltaT);
        writeFile << currentTime << "," << numS << "," << numE << "," << numA << "," << numI << "," << numR << "," << numQS << "," << numQE << "," << numQA << "," << numQI << "," << numQR << "," << numCR << "\n";
    }
    writeFile.close();
    const double ratioCR = numCR / (double)networkSize;

    if (deletion && ratioCR < 0.01) {
        std::cout << "Deleting file with CR=" << ratioCR << "\n";
        CSV::deleteFile(writeFileName);
    }
    return ratioCR;
}  //* End of function KM::syncRun

//* Update one step for single node
void asyncUpdate() {
    //* Calculate total transition rate and Delta t
    double totalTransitionRate = 0.0;
    for (const unsigned& index : reactingIndex) {
        totalTransitionRate += nodes[index].m_transitionRate;
    }
    const double deltaT = std::log(1.0 / probabilityDistribution(randomEngine)) / totalTransitionRate;
    totalTransitionRate *= probabilityDistribution(randomEngine);

    //* Update Time
    currentTime += deltaT;

    //* Release from quarantine
    release(deltaT);

    //* Choose target node to be reacted
    std::vector<unsigned> shuffledReactingIndex(reactingIndex.begin(), reactingIndex.end());
    std::shuffle(shuffledReactingIndex.begin(), shuffledReactingIndex.end(), randomEngine);
    double cumulativeTransitionRate = 0.0;
    unsigned target;
    for (unsigned i = 0; i < shuffledReactingIndex.size(); ++i) {
        target = shuffledReactingIndex[i];
        cumulativeTransitionRate += nodes[target].m_transitionRate;
        if (cumulativeTransitionRate > totalTransitionRate) {
            break;
        }
    }

    //* React target node and update transition rate
    const int intState = stateToInt.at(nodes[target].m_state);
    const double transitionRate = nodes[target].m_transitionRate;
    switch (intState) {
        //* S process
        case 0: {
            const double p = probabilityDistribution(randomEngine);
            //* S -> QS
            if (p * transitionRate < nodes[target].m_quarantineNeighbor * X_QX) {
                nodes[target].m_state = "QS";
                --numS;
                ++numQS;
                nodes[target].m_transitionRate = 0.0;
                reactingIndex.erase(target);
                nodes[target].m_quarantineTime = 0.0;
                quarantineIndex.emplace(target);
            }
            //* S -> E
            else {
                nodes[target].m_state = "E";
                --numS;
                ++numE;
                nodes[target].m_transitionRate = E_A + E_I + nodes[target].m_quarantineNeighbor * X_QX;
            }
            break;
        }
        //* E process
        case 1: {
            const double p = probabilityDistribution(randomEngine);
            //* E -> A
            if (p * transitionRate < E_A) {
                nodes[target].m_state = "A";
                --numE;
                ++numA;
                nodes[target].m_transitionRate = A_R + nodes[target].m_quarantineNeighbor * X_QX;
                for (const unsigned& neighbor : nodes[target].m_neighbors) {
                    if (nodes[neighbor].m_state == "S") {
                        nodes[neighbor].m_transitionRate += S_E;
                        reactingIndex.emplace(neighbor);
                    }
                }
            }
            //* E -> I
            else if (p * transitionRate < E_A + E_I) {
                nodes[target].m_state = "I";
                --numE;
                ++numI;
                nodes[target].m_transitionRate = I_QI + nodes[target].m_quarantineNeighbor * X_QX;
                for (const unsigned& neighbor : nodes[target].m_neighbors) {
                    if (nodes[neighbor].m_state == "S") {
                        nodes[neighbor].m_transitionRate += S_E;
                        reactingIndex.emplace(neighbor);
                    }
                }
            }
            //* E -> QE
            else {
                nodes[target].m_state = "QE";
                --numE;
                ++numQE;
                nodes[target].m_transitionRate = E_A + E_I;
                nodes[target].m_quarantineTime = 0.0;
                quarantineIndex.emplace(target);
            }
            break;
        }
        //* A process
        case 2: {
            const double p = probabilityDistribution(randomEngine);
            //* A -> R
            if (p * transitionRate < A_R) {
                nodes[target].m_state = "R";
                --numA;
                ++numR;
                nodes[target].m_transitionRate = X_QX * nodes[target].m_quarantineNeighbor;
                if (nodes[target].m_transitionRate < 1e-5) {
                    nodes[target].m_transitionRate = 0.0;
                    reactingIndex.erase(target);
                }
                for (const unsigned& neighbor : nodes[target].m_neighbors) {
                    if (nodes[neighbor].m_state == "S") {
                        nodes[neighbor].m_transitionRate -= S_E;
                        if (nodes[neighbor].m_transitionRate < 1e-5) {
                            nodes[neighbor].m_transitionRate = 0.0;
                            reactingIndex.erase(neighbor);
                        }
                    }
                }
            }
            //* A -> QA
            else {
                nodes[target].m_state = "QA";
                --numA;
                ++numQA;
                nodes[target].m_transitionRate = A_R;
                for (const unsigned& neighbor : nodes[target].m_neighbors) {
                    ++nodes[neighbor].m_quarantineNeighbor;
                    if (nodes[neighbor].m_state == "S") {
                        nodes[neighbor].m_transitionRate += X_QX - S_E;
                    } else if (nodes[neighbor].m_state == "R") {
                        reactingIndex.emplace(neighbor);
                        nodes[neighbor].m_transitionRate += X_QX;
                    } else if (nodes[neighbor].m_state == "E" || nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                        nodes[neighbor].m_transitionRate += X_QX;
                    }
                }
            }
            break;
        }
        //* I process
        case 3: {
            nodes[target].m_state = "QI";
            --numI;
            ++numQI;
            nodes[target].m_transitionRate = QI_CR;
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                ++nodes[neighbor].m_quarantineNeighbor;
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate += X_QX - S_E;
                } else if (nodes[neighbor].m_state == "R") {
                    reactingIndex.emplace(neighbor);
                    nodes[neighbor].m_transitionRate += X_QX;
                } else if (nodes[neighbor].m_state == "E" || nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                    nodes[neighbor].m_transitionRate += X_QX;
                }
            }
            break;
        }
        //* R process
        case 4: {
            nodes[target].m_state = "QR";
            --numR;
            ++numQR;
            nodes[target].m_transitionRate = 0.0;
            reactingIndex.erase(target);
            nodes[target].m_quarantineTime = 0.0;
            quarantineIndex.emplace(target);
            break;
        }
        //* QE process
        case 6: {
            const double p = probabilityDistribution(randomEngine);
            //* QE -> QA
            if (p * transitionRate < E_A) {
                nodes[target].m_state = "QA";
                --numQE;
                ++numQA;
                nodes[target].m_transitionRate = A_R;
                for (const unsigned& neighbor : nodes[target].m_neighbors) {
                    ++nodes[neighbor].m_quarantineNeighbor;
                    if (nodes[neighbor].m_state == "S" || nodes[neighbor].m_state == "R") {
                        nodes[neighbor].m_transitionRate += X_QX;
                        reactingIndex.emplace(neighbor);
                    } else if (nodes[neighbor].m_state == "E" || nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                        nodes[neighbor].m_transitionRate += X_QX;
                    }
                }
            }
            //* QE -> QI
            else {
                nodes[target].m_state = "QI";
                --numQE;
                ++numQI;
                nodes[target].m_transitionRate = QI_CR;
                for (const unsigned& neighbor : nodes[target].m_neighbors) {
                    ++nodes[neighbor].m_quarantineNeighbor;
                    if (nodes[neighbor].m_state == "S" || nodes[neighbor].m_state == "R") {
                        nodes[neighbor].m_transitionRate += X_QX;
                        reactingIndex.emplace(neighbor);
                    } else if (nodes[neighbor].m_state == "E" || nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                        nodes[neighbor].m_transitionRate += X_QX;
                    }
                }
            }
            break;
        }
        //* QA process
        case 7: {
            nodes[target].m_state = "CR";
            --numQA;
            ++numCR;
            nodes[target].m_transitionRate = 0.0;
            reactingIndex.erase(target);
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                --nodes[neighbor].m_quarantineNeighbor;
                if (nodes[neighbor].m_state == "S" || nodes[neighbor].m_state == "R") {
                    nodes[neighbor].m_transitionRate -= X_QX;
                    if (nodes[neighbor].m_transitionRate < 1e-5) {
                        reactingIndex.erase(neighbor);
                    }
                } else if (nodes[neighbor].m_state == "E" || nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                    nodes[neighbor].m_transitionRate -= X_QX;
                }
            }
            break;
        }
        //* QI Process
        case 8: {
            nodes[target].m_state = "CR";
            --numQI;
            ++numCR;
            nodes[target].m_transitionRate = 0.0;
            reactingIndex.erase(target);
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                --nodes[neighbor].m_quarantineNeighbor;
                if (nodes[neighbor].m_state == "S" || nodes[neighbor].m_state == "R") {
                    nodes[neighbor].m_transitionRate -= X_QX;
                    if (nodes[neighbor].m_transitionRate < 1e-5) {
                        reactingIndex.erase(neighbor);
                    }
                } else if (nodes[neighbor].m_state == "E" || nodes[neighbor].m_state == "A" || nodes[neighbor].m_state == "I") {
                    nodes[neighbor].m_transitionRate -= X_QX;
                }
            }
            break;
        }
        //* QS, QR, CR process
        default: {
            break;
        }
    }
}  //* End of function KM::asyncUpdate

double asyncRun() {
    //* Define write file
    const std::string writeFileName = directory + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS << "," << numE << "," << numA << "," << numI << "," << numR << "," << numQS << "," << numQE << "," << numQA << "," << numQI << "," << numQR << "," << numCR << "\n";
    while (!reactingIndex.empty()) {
        asyncUpdate();
        writeFile << currentTime << "," << numS << "," << numE << "," << numA << "," << numI << "," << numR << "," << numQS << "," << numQE << "," << numQA << "," << numQI << "," << numQR << "," << numCR << "\n";
    }

    writeFile.close();
    const double ratioCR = numCR / (double)networkSize;

    if (deletion && ratioCR < 0.01) {
        std::cout << "Deleting file with CR=" << ratioCR << "\n";
        CSV::deleteFile(writeFileName);
    }
    return ratioCR;
}  //* End of function KM::asyncRun
}  // namespace KM
