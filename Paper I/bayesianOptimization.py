from importedLibraries import *

'''CODE FOR EXPECTED IMPROVEMENT METHOD'''

class GetNextIndexExpectedImprovement(object):
    def __init__(self, importedData, predictedDataObj, indexesToCalculate, thresholdEstimate=False):
        self.importedData = importedData
        self.y_gpr = predictedDataObj.y_gpr
        self.y_std = predictedDataObj.y_std
        self.energyScaler = predictedDataObj.energyScaler
        self.indexesToCalculate = indexesToCalculate
        self.alpha = 0.05
        self.setYforSelectingIndex()
        self.thresholdEstimateGuess = max(self.yForSelectingIndex) if not thresholdEstimate else thresholdEstimate
        #If thresholdEstimate has a poor initial condition, the algorithm gets stuck
        if thresholdEstimate > max(self.yForSelectingIndex) or thresholdEstimate < 0:
            self.thresholdEstimateGuess = max(self.yForSelectingIndex)
        self.thresholdEstimate = self.getThresholdEstimate()
        self.nextIndex = self.getNextIndexToCalculateMinimization()
        self.FNestimate = self.getFNestimate(self.thresholdEstimate)
        self.TPestimate = self.getTPestimate(self.thresholdEstimate)
        self.FNrate = self.getFalseNegativeRate()

    def setYforSelectingIndex(self):
        self.yForSelectingIndex = -self.y_gpr

    def getFNestimate(self, thresholdEstimate):
        FNestimate = 0
        for i in range(0, len(self.yForSelectingIndex)):
            if i not in self.indexesToCalculate:
                mean, std = self.yForSelectingIndex[i], self.y_std[i]
                if std > 0:
                    Z = (mean - thresholdEstimate)/std
                    #print(Z)
                    FNestimate = FNestimate + norm.cdf(Z)
                    #print(norm.cdf(Z))
        return FNestimate

    def getTPestimate(self, thresholdEstimate):
        TPestimate = 0
        for i in range(0, len(self.yForSelectingIndex)):
            if i in self.indexesToCalculate and self.yForSelectingIndex[i] > thresholdEstimate:
                TPestimate = TPestimate + 1
        return TPestimate

    def getFalseNegativeRate(self):
        return self.FNestimate/(self.FNestimate + self.TPestimate)

    def thresholdEstimateMinimizationFunction(self, thresholdEstimate):
        FNestimate = self.getFNestimate(thresholdEstimate)
        TPestimate = self.getTPestimate(thresholdEstimate)
        self.minimizationFunction = FNestimate + TPestimate - len(self.yForSelectingIndex)*self.alpha
        return self.minimizationFunction

    def getThresholdEstimate(self):
        return fsolve(self.thresholdEstimateMinimizationFunction, self.thresholdEstimateGuess)

    def getNextIndexToCalculate(self):
        acqFunctionList = []
        for i in range(len(self.yForSelectingIndex)):
            mean, std = self.yForSelectingIndex[i], self.y_std[i]
            if std > 0 and i not in self.indexesToCalculate:
                Z = (mean - self.thresholdEstimate)/std
                acquisitionFunction = (mean - self.thresholdEstimate)*norm.cdf(Z) + std*norm.pdf(Z)
                acqFunctionList.append(acquisitionFunction)
            else:
                acqFunctionList.append(-1e10)
        nextIndex = acqFunctionList.index(max(acqFunctionList))
        if nextIndex in self.indexesToCalculate:
            nextIndex = self.getRandomNewIndex()
        return nextIndex

    def getNextIndexToCalculateMinimization(self):
        acqFunctionList = []
        minMean = min(self.y_gpr)
        for i in range(len(self.y_gpr)):
            mean, std = self.y_gpr[i], self.y_std[i]
            if std > 0 and i not in self.indexesToCalculate:
                Z = (mean - minMean)/std
                acquisitionFunction = (mean - minMean)*(1 - norm.cdf(Z)) - std*(norm.pdf(Z))
                acqFunctionList.append(acquisitionFunction)
            else:
                acqFunctionList.append(1e10)
        nextIndex = acqFunctionList.index(min(acqFunctionList))
        self.expectedImprovement = abs(min(acqFunctionList))
        # print("self.expectedImprovement", self.expectedImprovement)
        if nextIndex in self.indexesToCalculate:
            nextIndex = self.getRandomNewIndex()
        return nextIndex

    def getRandomNewIndex(self):
        randomIdx = random.randint(0, len(self.importedData.DFTdataDict["xDirect"])-1)
        while randomIdx in self.indexesToCalculate:
            randomIdx = random.randint(0, len(self.importedData.DFTdataDict["xDirect"])-1)
        return randomIdx
