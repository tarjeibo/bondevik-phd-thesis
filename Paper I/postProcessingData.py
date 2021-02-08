from importedLibraries import *


class GatherDataExpectedImprovement(object):
    def __init__(self, importedData, predictedDataObj, nextIndexCalcObj, gatheredData=False, printData=False):
        self.importedData = importedData
        self.predictedDataObj = predictedDataObj
        self.nextIndexCalcObj = nextIndexCalcObj
        self.gatheredData = gatheredData
        self.appendData()
        if printData==True: self.printBayesianOptimizationVariables()

    def appendData(self):
        if self.gatheredData==False: #First iteration
            self.iterationCounterList = [float(len(self.nextIndexCalcObj.indexesToCalculate))]
            self.FNestimateList = [float(self.nextIndexCalcObj.FNestimate)]
            self.TPestimateList = [float(self.nextIndexCalcObj.TPestimate)]
            self.FNrateList = [float(self.nextIndexCalcObj.FNrate)]
            self.thresholdEstimateList = [float(self.nextIndexCalcObj.thresholdEstimate)]
            self.minimizationFunctionList = [float(self.nextIndexCalcObj.minimizationFunction)]
            self.y_gprList = [self.predictedDataObj.energyScaler.inverse_transform(self.predictedDataObj.y_gpr)]
            self.y_stdList = [self.nextIndexCalcObj.y_std]
            self.calculatedXDirectList = [self.importedData.DFTdataDict["xDirect"][self.nextIndexCalcObj.indexesToCalculate]]
            self.calculatedEnergyList = [self.nextIndexCalcObj.energyScaler.inverse_transform(self.nextIndexCalcObj.y_gpr)[self.nextIndexCalcObj.indexesToCalculate]]
            self.expectedImprovementList = [self.nextIndexCalcObj.expectedImprovement]
            self.minYgprIndexList = [np.argmin(self.nextIndexCalcObj.y_gpr)]
        else: #Setting instance variable to previously gathered data, then gathering the new data point
            self.iterationCounterList = self.gatheredData.iterationCounterList
            self.iterationCounterList.append(float(len(self.nextIndexCalcObj.indexesToCalculate)))
            self.FNestimateList = self.gatheredData.FNestimateList
            self.FNestimateList.append(float(self.nextIndexCalcObj.FNestimate))
            self.TPestimateList = self.gatheredData.TPestimateList
            self.TPestimateList.append(float(self.nextIndexCalcObj.TPestimate))
            self.FNrateList = self.gatheredData.FNrateList
            self.FNrateList.append(float(self.nextIndexCalcObj.FNrate))
            self.thresholdEstimateList = self.gatheredData.thresholdEstimateList
            self.thresholdEstimateList.append(float(self.nextIndexCalcObj.thresholdEstimate))
            self.minimizationFunctionList = self.gatheredData.minimizationFunctionList
            self.minimizationFunctionList.append(float(self.nextIndexCalcObj.minimizationFunction))
            self.y_gprList = self.gatheredData.y_gprList
            self.y_gprList.append(self.predictedDataObj.energyScaler.inverse_transform(self.predictedDataObj.y_gpr))
            self.y_stdList = self.gatheredData.y_stdList
            self.y_stdList.append(self.nextIndexCalcObj.y_std)
            self.calculatedXDirectList = self.gatheredData.calculatedXDirectList
            self.calculatedXDirectList.append(self.importedData.DFTdataDict["xDirect"][self.nextIndexCalcObj.indexesToCalculate])
            self.calculatedEnergyList = self.gatheredData.calculatedEnergyList
            self.calculatedEnergyList.append(self.nextIndexCalcObj.energyScaler.inverse_transform(self.nextIndexCalcObj.y_gpr)[self.nextIndexCalcObj.indexesToCalculate])
            self.expectedImprovementList = self.gatheredData.expectedImprovementList
            self.expectedImprovementList.append(self.nextIndexCalcObj.expectedImprovement)
            self.minYgprIndexList = self.gatheredData.minYgprIndexList
            self.minYgprIndexList.append(np.argmin(self.nextIndexCalcObj.y_gpr))

    def printBayesianOptimizationVariables(self):
        print("FN:", "%.2f" % float(self.nextIndexCalcObj.FNestimate), "TP:", self.nextIndexCalcObj.TPestimate, "FNrate:", float(self.nextIndexCalcObj.FNrate), "ThresholdEstimate:", "%.2f" % float(self.nextIndexCalcObj.thresholdEstimate), "Max y_gpr:", "%.2f" % max(self.nextIndexCalcObj.y_gpr), "Min.func:", "%.6f" % float(self.nextIndexCalcObj.minimizationFunction))

    def writeFNrateDataToCsv(self, fileNameIterator):
        dataToCsv = [self.iterationCounterList, self.FNestimateList, self.TPestimateList, self.FNrateList, self.thresholdEstimateList, self.minimizationFunctionList]
        fileNameWithLocation = 'ExportedData/expectedImprovementFNrateData' + str(fileNameIterator) + '.csv'
        with open(fileNameWithLocation, 'w', newline='') as myFile:
            wr = csv.writer(myFile)
            wr.writerows(dataToCsv)

    def writeYDataToCsv(self, fileNameIterator):
        dataTypes = ["calculatedEnergy", "calculatedXDirect", "Ygpr", "Ystd"]
        numberOfCalculatedDataPts = len(self.calculatedXDirectList)
        for dataType in dataTypes:
            y_grpData = self.y_gprList
            if dataType == "calculatedEnergy":
                dataToCsv = self.calculatedEnergyList
            elif dataType == "calculatedXDirect":
                dataToCsv = self.calculatedXDirectList
            elif dataType == "Ygpr":
                dataToCsv = y_grpData
            else:
                dataToCsv = []
                for i in range(0, numberOfCalculatedDataPts):
                    dataToCsv.append(self.y_stdList[i]*np.std(y_grpData[i]))
            fileNameWithLocation = 'ExportedData/expectedImprovementPredicted' + str(dataType) + str(fileNameIterator) + '.csv'
            print("Saved the following file: ", fileNameWithLocation)
            with open(fileNameWithLocation, 'w', newline='') as myFile:
                wr = csv.writer(myFile)
                wr.writerows(dataToCsv)

    # def writeConvergenceDataToCsv(self, fileNameIterator, convergenceCriterionList=False):
    #     if convergenceCriterionList == False:
    #         convergenceCriterionList = [1e-10, 1e-20, 1e-30, 1e-40, 1e-50, 1e-60]
    #     numberOfCalculationsAtConvergenceList, minimumEnergyIndexAtConvergenceList = self.getSingleRunConvergenceData(convergenceCriterionList)
    #     dataToCsv = [convergenceCriterionList, numberOfCalculationsAtConvergenceList, minimumEnergyIndexAtConvergenceList]
    #     fileNameWithLocation = 'ExportedData/expectedImprovementConvergenceData' + str(fileNameIterator) + '.csv'
    #     print("Saved the following file: ", fileNameWithLocation)
    #     with open(fileNameWithLocation, 'w', newline='') as myFile:
    #         wr = csv.writer(myFile)
    #         wr.writerows(dataToCsv)

    def getCPUcost(self, i):
        initializationPoints = len(self.nextIndexCalcObj.indexesToCalculate) - len(self.iterationCounterList) + 1
        calculatedIndexes = self.nextIndexCalcObj.indexesToCalculate[:i + initializationPoints]
        totalCPUcost = sum(self.importedData.DFTdataDict["energyCPU"])
        CPUcostFromCalculatedIndexes = sum(self.importedData.DFTdataDict["energyCPU"][calculatedIndexes])
        CPUcostFromDescriptors = 0
        for descriptor in self.importedData.descriptors:
            if descriptor == "prePES" or descriptor == "nsw3" or descriptor == "nsw6":
                CPUcostFromDescriptors = CPUcostFromDescriptors + sum(self.importedData.DFTdataDict[str(descriptor + "CPU")])
        fractionalCPUcost = (CPUcostFromCalculatedIndexes + CPUcostFromDescriptors)/totalCPUcost
        return fractionalCPUcost

    def getSingleRunConvergenceData(self, convergenceCriterionList=False):
        if convergenceCriterionList == False:
            convergenceCriterionList = [1e-10, 1e-20, 1e-30, 1e-40, 1e-50, 1e-60]
        numberOfCalculationsAtConvergenceList = []
        minimumEnergyIndexAtConvergenceList = []
        CPUcostAtConvergenceList = []
        for convergenceCriterion in convergenceCriterionList:
            for i, expectedImprovement in enumerate(self.expectedImprovementList):
                if expectedImprovement < convergenceCriterion:
                    numberOfCalculationsAtConvergenceList.append(len(self.calculatedEnergyList[i]))
                    minimumEnergyIndexAtConvergenceList.append(self.minYgprIndexList[i])
                    CPUcostAtConvergenceList.append(self.getCPUcost(i))
                    break
                elif i == len(self.expectedImprovementList) - 1:
                    numberOfCalculationsAtConvergenceList.append(-1)
                    minimumEnergyIndexAtConvergenceList.append(-1)
                    CPUcostAtConvergenceList.append(-1)
        return convergenceCriterionList, numberOfCalculationsAtConvergenceList, minimumEnergyIndexAtConvergenceList, CPUcostAtConvergenceList

class ConvergenceData(object):
    def __init__(self, convergenceDataObj=False, convergenceCriterionList=False, singleRunConvergenceData=False):
        self.convergenceDataObj = convergenceDataObj
        self.convergenceCriterionList = convergenceCriterionList
        self.singleRunConvergenceData = singleRunConvergenceData
        self.gatherConvergenceData()

    def gatherConvergenceData(self):
        if self.convergenceDataObj == False: #First iteration
            self.convergenceData = []
        else:
            self.convergenceData = self.convergenceDataObj.convergenceData
            self.convergenceData.append(self.singleRunConvergenceData)

    def writeConvergenceDataToCsv(self, descriptorIdx=False):
        numberOfConvergenceCriterions = len(self.convergenceData[0][0])
        numberOfRuns = len(self.convergenceData)
        dataToCsv = []
        for i in range(0, numberOfConvergenceCriterions):
            numberOfCalculationsAtConvergenceList = [self.convergenceCriterionList[i]]
            minimumEnergyIndexAtConvergenceList = [self.convergenceCriterionList[i]]
            CPUcostAtConvergenceList = [self.convergenceCriterionList[i]]
            for j in range(0, numberOfRuns):
                numberOfCalculationsAtConvergenceList.append(self.convergenceData[j][1][i])
                minimumEnergyIndexAtConvergenceList.append(self.convergenceData[j][2][i])
                CPUcostAtConvergenceList.append(self.convergenceData[j][3][i])
            dataToCsv.append(numberOfCalculationsAtConvergenceList)
            dataToCsv.append(minimumEnergyIndexAtConvergenceList)
            dataToCsv.append(CPUcostAtConvergenceList)
        fileName = 'expectedImprovementConvergenceData' + '.csv'
        print("Saved the following file: ", fileName)
        with open(fileName, 'w', newline='') as myFile:
            wr = csv.writer(myFile)
            wr.writerows(dataToCsv)
