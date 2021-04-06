import pandas as pd
import numpy as np

rootDirectory = "../data/KM/"


def readCSV(t_fileName):
    data = pd.read_csv(t_fileName, sep=',', header=None)
    data = data.values.transpose()
    if (len(data) == 1):
        return data[0]
    else:
        return tuple([row for row in data])



#* Read Data
# def readKM(network_type, network_size, mean_degree, S_E, E_AI, pA, I_QI, A_R, QI_CR, X_QX, tau, t_randomEngineSeed=None):
#     directory = rootDirectory + network_type + ",N{:.1e},M{:d}/".format(network_size, mean_degree)
#     fileName = "SE{:.2f},EAI{:.2f},pA{:.2f},IQI{:.2f},AR{:.2f},QICR{:.2f},XQX{:.2f},T{:.2f}".format(S_E, E_AI, pA, I_QI, A_R, QI_CR, X_QX, tau)
#     if (t_randomEngineSeed == None):
#         fullFileName = directory + fileName + ".txt"
#     else:
#         fullFileName = directory + fileName + "-{:d}.txt".format(int(t_randomEngineSeed))

#     return readCSV(fullFileName)

def readKM(fileName):
    return readCSV("../data/KM/" + fileName)

def readReal():
    return readCSV("realConfirmed_208.txt")


def readFig4():
    return readCSV("fig4_t_S_QS_E_QE_A_QA_I_QI_R_QR_CR.txt")

if __name__ == "__main__":
    print("This is a moduel parameters.py")
