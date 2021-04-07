#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../library/stringFormat.hpp"

const std::string networkName(const std::string& t_networkType, const unsigned& t_networkSize, const unsigned& t_meanDegree) {
    return t_networkType + ",N" + to_stringWithExponent((double)t_networkSize, 1) + ",M" + std::to_string(t_meanDegree) + "/";
}

const std::string rateName(const std::vector<double>& t_rates, const int& t_randomEngineSeed) {
    if (t_rates.size() != 8) {
        std::cout << "In rate name function, input rate vector is not size of 8\n";
        exit(1);
    }

    const std::string fileName = "SE" + to_stringWithPrecision(t_rates[0], 3) + ",EAI" + to_stringWithPrecision(t_rates[1], 3) + ",pA" + to_stringWithPrecision(t_rates[2], 3) + ",IQI" + to_stringWithPrecision(t_rates[3], 3) + ",AR" + to_stringWithPrecision(t_rates[4], 3) + ",QICR" + to_stringWithPrecision(t_rates[5], 3) + ",XQX" + to_stringWithPrecision(t_rates[6], 3) + ",T" + to_stringWithPrecision(t_rates[7], 3);

    return t_randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_randomEngineSeed) + ".txt";
}

