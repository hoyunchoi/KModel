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
def KM(fileName):
    intToState = {0: "S", 1: "E", 2: "A", 3: "I", 4: "R", 5: "QS", 6: "QE", 7: "QA", 8: "QI", 9: "QR", 10: "CR"}
    df = pd.read_csv("../data/KM/" + fileName, header=None)
    data = {}
    for index, value in enumerate(df.values.transpose()):
        data[intToState[index]] = value
    return data


def KM_Rt(fileName):
    intToState = {0: "S", 1: "E", 2: "A", 3: "I", 4: "R", 5: "QS", 6: "QE", 7: "QA", 8: "QI", 9: "QR", 10: "CR", 11: "L", 12: "Rt"}
    df = pd.read_csv("../data/KM_Rt/" + fileName, header=None)
    data = {}
    for index, value in enumerate(df.values.transpose()):
        data[intToState[index]] = value
    return data


def real_confirmed(t=209):
    return readCSV("../data/COVID/realData/confirmed_{:d}.txt".format(t))


def real_Rt(t=342):
    return readCSV("../data/COVID/realData/Rt_{:d}.txt".format(t))


def fig4():
    intToState = {0: "t", 1: "S", 2: "QS", 3: "E", 4: "QE", 5: "A", 6: "QA", 7: "I", 8: "QI", 9: "R", 10: "QR", 11: "CR"}
    df = pd.read_csv("../data/COVID/realData/fig4_data.txt", header=None)
    data = {}
    for index, value in enumerate(df.values.transpose()):
        data[intToState[index]] = value[::100]
    return data


def getMinimumEnergy(networkSize):
    minimum_energy = 1e6
    minimum_line = ""
    with open("energyList.txt") as file:
        for line in file.readlines():
            if "{:.1e}".format(networkSize) in line:
                energy = float(line[line.find(": ") + 2:])
                if energy < minimum_energy:
                    minimum_energy = energy
                    minimum_line = line
    return minimum_line


if __name__ == "__main__":
    print("This is a moduel read_data.py")
