from importedLibraries import *


class PlotDataExpectedImprovement(object):
    def __init__(self, importedData, predictedDataObj, nextIndexCalcObj, gatheredData):
        self.importedData = importedData
        self.predictedDataObj = predictedDataObj
        self.nextIndexCalcObj = nextIndexCalcObj
        self.gatheredData = gatheredData
        self.xDirect = self.importedData.DFTdataDict["xDirect"]
        self.yDirect = self.importedData.DFTdataDict["yDirect"]


    def plot1dEnergyData(self, plotYsamples=False):
        plt.figure()

        y_gpr_unScaled = self.predictedDataObj.energyScaler.inverse_transform(self.predictedDataObj.y_gpr)
        y_std_unScaled = self.predictedDataObj.y_std*np.std(y_gpr_unScaled)

        ySelected = self.predictedDataObj.energyScaler.inverse_transform(self.predictedDataObj.y)
        Xselected = self.predictedDataObj.trainingDataScalerDict[0].inverse_transform(self.predictedDataObj.X[:,0])
        lw = 2
        plt.scatter(self.xDirect, self.importedData.DFTdataDict['unScaledEnergy'], c='k', label='Data points')

        plt.plot(self.xDirect, y_gpr_unScaled, color='darkorange', lw=lw, label='Gaussian Process Fit')
        plt.fill_between(self.xDirect, y_gpr_unScaled - y_std_unScaled, y_gpr_unScaled + y_std_unScaled, color='darkorange', alpha=0.2, label="Uncertainty")
        plt.scatter(self.xDirect[self.nextIndexCalcObj.indexesToCalculate], y_gpr_unScaled[self.nextIndexCalcObj.indexesToCalculate], c='r', s=100, label='Calculated data points')

        plt.ylabel('Energy / eV', fontsize=14)
        plt.xlabel("$\mathit{x}$", fontsize=18)
        plt.legend(loc="best",  scatterpoints=1, prop={'size': 10})













#some more space
