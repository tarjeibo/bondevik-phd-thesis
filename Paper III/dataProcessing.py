#Module to that plots and saves data (and maybe imports, as well)

from physicalConstants import *
from importedLibraries import *

class ProcessData(object):
    def __init__(self, grainBoundaryObj, gridParamsObj, chargeDensityObj, poissonSolverObj):
        self.grainBoundaryObj = grainBoundaryObj
        self.gridParamsObj = gridParamsObj
        self.chargeDensityObj = chargeDensityObj
        self.poissonSolverObj = poissonSolverObj
        self.solverType = self.poissonSolverObj.solverType
        self.setXandConcentrationAndPotential()


    def setXandConcentrationAndPotential(self):
        halfGBwidth = self.gridParamsObj.systemWidth/2
        self.setXranges(halfGBwidth)
        if self.solverType == "Jacobi":
            SCLconcentrationsJacobi = self.poissonSolverObj.chargeDensityObj.getSCLconcentrationsJacobi()
            self.c_V_list = np.array(SCLconcentrationsJacobi.c_V_SCL_list)
            self.c_OH_list = np.array(SCLconcentrationsJacobi.c_OH_SCL_list)
            self.c_Acc_list = np.array(SCLconcentrationsJacobi.c_Acc_SCL_list)
            self.c_AccV_list = np.array(SCLconcentrationsJacobi.c_AccV_SCL_list)
            self.c_AccOH_list = np.array(SCLconcentrationsJacobi.c_AccOH_SCL_list)
            self.phiRange = self.poissonSolverObj.phiList
        elif self.solverType == "BVP":
            self.c_V_list = self.reverseAndMergeArray(np.array(self.poissonSolverObj.mainBVPresults.c_V_SCLlist))
            self.c_OH_list = self.reverseAndMergeArray(np.array(self.poissonSolverObj.mainBVPresults.c_OH_SCLlist))
            self.c_Acc_list = self.reverseAndMergeArray(np.array(self.poissonSolverObj.mainBVPresults.c_Acc_SCLlist))
            self.convertConcentrationListsToOneStepPerAtomicLayer()
            self.c_V_list[len(self.c_V_list)//2] = self.poissonSolverObj.chargeDensityObj.coreConcentrations.c_V_core
            self.c_OH_list[len(self.c_OH_list)//2] = self.poissonSolverObj.chargeDensityObj.coreConcentrations.c_OH_core
            self.phiRange = self.reverseAndMergeArray(np.array(self.poissonSolverObj.mainBVPresults.phiSCLlist))
        elif self.solverType == "Analytical":
            analyticalSCLlength = self.poissonSolverObj.mainAnalyticalResults.SCLlength
            self.xRangePotential = np.linspace(-analyticalSCLlength, analyticalSCLlength, self.gridParamsObj.numberOfSteps)
            self.phiRange = self.reverseAndMergeArray(np.array(self.poissonSolverObj.mainAnalyticalResults.phiSCLlist))

    def setXranges(self, halfGBwidth):
        self.xRangeConcentration = np.linspace(-halfGBwidth, halfGBwidth, self.gridParamsObj.numberOfLayers)
        self.xRangePotential = np.linspace(-halfGBwidth, halfGBwidth, self.gridParamsObj.numberOfSteps)

    def reverseAndMergeArray(self, array):
        tempArray1 = array
        tempArray2 = np.delete(tempArray1, 0)
        reversedArray = np.flipud(tempArray2)
        return np.append(reversedArray, array)

    def convertConcentrationListsToOneStepPerAtomicLayer(self):
        c_V_list_temp = self.c_V_list
        c_OH_list_temp = self.c_OH_list
        c_Acc_list_temp = self.c_Acc_list
        self.c_V_list, self.c_OH_list, self.c_Acc_list = [], [], []
        for i in range(0, len(self.xRangeConcentration)):
            j = (np.abs(self.xRangePotential - self.xRangeConcentration[i]).argmin())
            self.c_V_list.append(c_V_list_temp[j])
            self.c_OH_list.append(c_OH_list_temp[j])
            self.c_Acc_list.append(c_Acc_list_temp[j])

    def plotConcentrations(self, savePlot=False):
        label = self.poissonSolverObj.solverType
        if label == "Analytical":
            return
        plt.plot(self.xRangeConcentration, self.c_V_list, label="Vacancies " + str(label))
        plt.plot(self.xRangeConcentration, self.c_OH_list, label="Protons " + str(label))
        if label == "Jacobi":
            plt.plot(self.xRangeConcentration, self.c_Acc_list, label="Acc-doping " + str(label))
            plt.plot(self.xRangeConcentration, self.c_AccV_list, label="AccV-cluster " + str(label))
            plt.plot(self.xRangeConcentration, self.c_AccOH_list, label="AccOH-cluster " + str(label))
        plt.legend(bbox_to_anchor=(1, 0.5), loc='upper left').get_frame().set_linewidth(0)
        plt.xlabel("x / m")
        plt.ylabel("Log concentrations / m$^{-3}$")
        plt.yscale("log")
        if savePlot == True:
            plt.savefig('fooChg.png', dpi=300, bbox_inches='tight')
        #plt.show()
        #plt.close()

    def plotPotential(self, savePlot=False):
        label = self.poissonSolverObj.solverType
        plt.plot(self.xRangePotential, self.phiRange, label=str(label))
        plt.legend(bbox_to_anchor=(1, 0.5), loc='upper left').get_frame().set_linewidth(0)
        plt.xlabel("x / m")
        plt.ylabel("Potential / V")
        if savePlot == True:
            plt.savefig('fooPot.png', dpi=300, bbox_inches='tight')
        #plt.show()
        #plt.close()

    def writeDataToCsv(self, fileName="PoissonCsvData"):
        with open(fileName, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            if self.solverType == "Jacobi":
                wr.writerows([self.xRangeConcentration, self.c_V_list, self.c_OH_list, self.c_Acc_list, self.c_AccV_list, self.c_AccOH_list, self.xRangePotential, self.phiRange])
            elif self.solverType == "BVP":
                wr.writerows([self.xRangeConcentration, self.c_V_list, self.c_OH_list, self.xRangePotential, self.phiRange])
            elif self.solverType == "Analytical":
                wr.writerows([self.xRangePotential, self.phiRange])


class ElectricalProperties(object):
    """docstring for ElectricalProperties."""
    def __init__(self, processedDataObj):
        self.processedDataObj = processedDataObj
        self.xRangeConcentration = self.processedDataObj.xRangeConcentration
        self.xRangePotential = self.processedDataObj.xRangePotential
        self.c_V_list = self.processedDataObj.c_V_list
        self.c_OH_list = self.processedDataObj.c_OH_list
        self.c_V_bulk = self.processedDataObj.grainBoundaryObj.bulkConcentrations.c_V_bulk
        self.c_OH_bulk = self.processedDataObj.grainBoundaryObj.bulkConcentrations.c_OH_bulk
        self.phiRange = self.processedDataObj.phiRange
        self.wRange = np.logspace(-5, 6, 4000)
        self.RgbOverRbulkRatio = self.getRgbOverRbulkRatio()
        self.solverType = self.processedDataObj.poissonSolverObj.solverType

    #Function assuming only protons, constant mobility and that segregated protons do not contribute to conductivity
    def getResistanceFromPotential(self):
        xLength = self.xRangeConcentration[-1] - self.xRangeConcentration[0]
        xStep = self.xRangePotential[1] - self.xRangePotential[0]
        self.RgbOverRbulkRatioFromPotentialList = []
        for i in range(0, len(self.xRangePotential)):
            RgbOverRbulkRatioAtX = np.exp((e*self.phiRange[i])/(kB*self.processedDataObj.grainBoundaryObj.T))*(xStep/xLength)
            self.RgbOverRbulkRatioFromPotentialList.append(RgbOverRbulkRatioAtX)
        return self.RgbOverRbulkRatioFromPotentialList

    #Function assumes system with only protons, and constant mobility through GB region
    def getRgbOverRbulkRatio(self):
        xLength = self.xRangeConcentration[-1] - self.xRangeConcentration[0]
        xStep = self.xRangeConcentration[1] - self.xRangeConcentration[0]
        self.RgbOverRbulkRatioList = []
        for i in range(0, len(self.xRangeConcentration)):
            RgbOverRbulkRatioAtX = (self.c_OH_bulk/xLength)*(xStep/self.c_OH_list[i])
            self.RgbOverRbulkRatioList.append(RgbOverRbulkRatioAtX)
        RgbOverRbulkRatio = sum(self.RgbOverRbulkRatioList)
        return RgbOverRbulkRatio

    #Function assumes system with only protons, and constant mobility through GB region
    #The capacitance is set to unity, the resistance to R_gb/R_bulk
    def getSingleGBimpedance(self):
        Z = 0
        j = sqrt(-1)
        RlistPotential = self.getResistanceFromPotential()
        for i in range(0, len(self.xRangePotential)):
            #R = self.RgbOverRbulkRatioList[i]
            R = RlistPotential[i]
            Z = Z + R/(1+R*self.wRange*j)
        return Z

    def plotSingleGBrealImpedance(self, savePlot=False):
        Z = self.getSingleGBimpedance()
        label = self.solverType
        plt.plot(self.wRange, Z.real, label=label)
        plt.xscale("log")
        plt.xlabel("w / 1/s")
        plt.ylabel("Z' (arbitrary units)")
        plt.legend(loc="best")
        if savePlot == True:
            plt.savefig('fooImp.png', dpi=300, bbox_inches='tight')

    def plotSingleGBimagImpedance(self, savePlot=False):
        Z = self.getSingleGBimpedance()
        ZimagNegative = [-x for x in Z.imag]
        label = self.solverType
        plt.plot(self.wRange, ZimagNegative, label=label)
        plt.xscale("log")
        plt.xlabel("w / 1/s")
        plt.ylabel("-Z'' (arbitrary units)")
        plt.legend(loc="best")
        if savePlot == True:
            plt.savefig('fooImp.png', dpi=300, bbox_inches='tight')


class SaveData(object):
    def __init__(self, electricalPropertiesObj):
        self.electricalPropertiesObj = electricalPropertiesObj
        self.saveSingleGBimpedanceToCsv()

    def saveSingleGBimpedanceToCsv(self):
        Z = self.electricalPropertiesObj.getSingleGBimpedance()
        dataToCsv = [self.electricalPropertiesObj.wRange, Z.real, Z.imag]
        fileName = 'singleGBimpedance' + str(self.electricalPropertiesObj.solverType) + '.csv'
        print("Saved the following file: ", fileName)
        with open(fileName, 'w', newline='') as myFile:
            wr = csv.writer(myFile)
            wr.writerows(dataToCsv)



def getXrowAndPhiRow(fileName, xRow, phiRow):
    with open(fileName, newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i == xRow:
                xRow = row
            if i == phiRow:
                phiRow = row
    return np.asarray(xRow).astype(float), np.asarray(phiRow).astype(float)

def calcConductivityRatio(fileName, xRow, phiRow, T, maxXrow):
    xRow, phiRow = getXrowAndPhiRow(fileName, xRow, phiRow)
    print("xRow", xRow)
    print("phiRow", phiRow)
    dX = xRow[1] - xRow[0]
    extraTermForModelM1 = 2*abs((max(xRow) - maxXrow))
    resistivityRatioSum = dX*np.sum(np.exp(e*phiRow/(kB*T))) + extraTermForModelM1
    sigmaGBoverSigmaBulk = 1/(resistivityRatioSum/(2*maxXrow))
    # print(sigmaGBoverSigmaBulk)
    # print("Did you remember to set the temperature correctly????")
    return sigmaGBoverSigmaBulk

def calc_conductivity_ratio_M3(fileName):
    xRow, c_OH_layer = getXrowAndPhiRow(fileName, 0, 2)
    # print("xRow", np.asarray(xRow))
    # print("c_OH_layer", np.asarray(c_OH_layer))
    c_OH_bulk = c_OH_layer[0]
    # print("c_OH_bulk", c_OH_bulk)
    return 1/np.average(c_OH_bulk/c_OH_layer)


if __name__ == '__main__':
    xRow = 0
    phiRow = 1
    # fileName = "T=600210GB,p=2.5e-05,T=600Analytical"
    fileName = "T=600cluster111GB"

    maxXrow = 4.890150934808013e-09 #111 GB, 100e-10 total length
    # maxXrow = 4.924263014047538e-09 #210 GB, 100e-10 total length

    T = 1200
    xRow, phiRow = 6, 7
    # print("sigma_ratio_old", calcConductivityRatio(fileName, xRow, phiRow, T, maxXrow))
    print('\n')
    print("sigma_ratio_new", calc_conductivity_ratio_M3(fileName))
