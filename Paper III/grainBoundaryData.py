#Module containing grain boundary data
from importedLibraries import *
from physicalConstants import *

class GrainBoundary(object):
    def __init__(self, hkl, a0, GBtype, segregationEnergies, approximateSystemWidth, T, dH, dS, pH2O, accDoping, dE_AccV=None, dE_AccOH=None, N_conf_cluster=None):
        self.hkl = hkl
        self.a0 = a0
        self.GBtype = GBtype
        self.segregationEnergies = segregationEnergies
        self.approximateSystemWidth = approximateSystemWidth
        self.T = T
        self.dH = dH
        self.dS = dS
        self.pH2O = pH2O
        self.accDoping = accDoping
        self.dE_AccV = dE_AccV
        self.dE_AccOH = dE_AccOH
        self.N_conf_cluster = N_conf_cluster
        self.numberOfLayers = self.getNumberOfLayers()
        self.layerLength = self.getLayerLength()
        self.layerArea = self.getLayerArea()
        self.unitCellVolume = self.getUnitCellVolume()
        self.atomicPlaneVolume = self.getAtomicPlaneVolume()
        print("self.layerArea", self.layerArea)
        print("self.atomicPlaneVolume", self.atomicPlaneVolume)
        self.setClusterEquilibriumConstants()
        self.setKappa()
        self.setBulkConcentrations()
        self.setSegregationEnergyLists()

    def setClusterEquilibriumConstants(self):
        if self.dE_AccV == None:
            self.K_AccV = 0
            self.K_AccOH = 0
        else:
            self.K_AccV = self.N_conf_cluster*np.exp(-self.dE_AccV/(kBev*self.T))
            self.K_AccOH = self.N_conf_cluster*np.exp(-self.dE_AccOH/(kBev*self.T))
        print("self.K_AccV, self.K_AccOH", self.K_AccV, self.K_AccOH)

    def setKappa(self):
        self.kappa = self.pH2O*np.exp((-self.dH/(kBev*self.T) + self.dS*1e-3/kBev))

    def getUnitCellVolume(self):
        return self.a0**3

    def getIntegerHkl(self):
        return int(self.hkl[0]), int(self.hkl[1]), int(self.hkl[2])

    def getLayerLength(self):
        h, k, l = self.getIntegerHkl()
        if self.hkl == '111': #Every other layer does not have any oxygen atoms
            return self.a0/(np.sqrt(h**2 + k**2 + l**2))
        elif self.hkl == '210':
            return self.a0/(2*np.sqrt(h**2 + k**2 + l**2))

    def getAtomicPlaneVolume(self):
        h, k, l = self.getIntegerHkl()
        if self.hkl == '111':
            return (np.sqrt(3)/6)*self.a0*np.sqrt(2)*self.a0**2 #length x area
        elif self.hkl == '210':
            return (self.a0**3)/2

    def getLayerArea(self):
        h, k, l = self.getIntegerHkl()
        if self.hkl == '111':
            return np.sqrt(2)*self.a0**2
        elif self.hkl == '210':
            return np.sqrt(5)*self.a0**2

    def getNumberOfLayers(self):
        return int(np.ceil((self.approximateSystemWidth/self.getLayerLength())//2)*2 + 1)

    def setSegregationEnergyLists(self):
        dE_V_list = self.createSegregationEnergyList(self.segregationEnergies.vacancySegregationEnergies)
        dE_OH_list = self.createSegregationEnergyList(self.segregationEnergies.protonSegregationEnergies)
        dE_AccV_list = self.createSegregationEnergyList(self.segregationEnergies.acceptorVacancySegregationEnergies)
        dE_AccOH_list = self.createSegregationEnergyList(self.segregationEnergies.acceptorProtonSegregationEnergies)
        self.segregationEnergyLists = SegregationEnergyLists(dE_V_list, dE_OH_list, dE_AccV_list, dE_AccOH_list)

    def createSegregationEnergyList(self, segregationEnergies):
        segregationEnergyList = []
        if self.hkl=='210':
            for i in range(0, self.getNumberOfLayers()):
                if i%2 == 0:
                    segregationEnergyList.append([0, 0])
                else:
                    segregationEnergyList.append([0])
        if self.hkl=='111':
            for i in range(0, self.getNumberOfLayers()):
                segregationEnergyList.append([0, 0, 0])
        if self.GBtype=="OneSidedInput":
            GBcoreIndex = self.getNumberOfLayers()//2
            for i, segregationEnergy in enumerate(segregationEnergies):
                segregationEnergyList[GBcoreIndex + i] = segregationEnergies[i]
                segregationEnergyList[GBcoreIndex - i] = segregationEnergies[i]
        elif self.GBtype=="TwoSidedInput":
            startIndex = self.getNumberOfLayers()//2 - len(segregationEnergies)//2
            for i, segregationEnergy in enumerate(segregationEnergies):
                segregationEnergyList[startIndex + i] = segregationEnergies[i]
        return segregationEnergyList

    def setBulkConcentrations(self):
        c_V_bulk, c_OH_bulk =  fsolve(self.bulkConcentrationsNumericalSolver, (0, 0))
        c_O_bulk = (c_OH_bulk**2)/(self.kappa*c_V_bulk)

        # print("CONCENTRATIONS ARE CHANGED MANUALLY INSIDE CODE!!!")
        # print("CONCENTRATIONS ARE CHANGED MANUALLY INSIDE CODE!!!")
        # print("CONCENTRATIONS ARE CHANGED MANUALLY INSIDE CODE!!!")
        # c_V_bulk, c_OH_bulk, c_O_bulk = 0, 0.1, 2.9

        c_AccV_bulk = -(c_V_bulk*self.K_AccV*(c_OH_bulk + 2*c_V_bulk))/(c_V_bulk*self.K_AccV - 1)
        c_Acc_bulk = 2*c_V_bulk + c_OH_bulk + c_AccV_bulk
        c_AccOH_bulk = self.accDoping - c_Acc_bulk - c_AccV_bulk
        self.bulkConcentrations = BulkConcentrations(c_V_bulk, c_OH_bulk, c_O_bulk, c_Acc_bulk, c_AccV_bulk, c_AccOH_bulk)
        print("c_V_bulk", c_V_bulk)
        print("c_OH_bulk", c_OH_bulk)
        print("Sum, O sites: ", c_V_bulk + c_OH_bulk + c_O_bulk + c_AccV_bulk + c_AccOH_bulk)
        print("Sum, Acc sites: ", c_Acc_bulk + c_AccV_bulk + c_AccOH_bulk)

    def bulkConcentrationsNumericalSolver(self, variables):
        c_V, c_OH = variables
        eqOne = 2*c_V + c_OH - self.accDoping*(1 - c_V*self.K_AccV)/(1 + self.K_AccV*c_V + self.K_AccOH*c_OH)
        eqTwo = c_OH**2 + (c_V*self.kappa)*(-3 + c_V + c_OH + (self.accDoping/(1 + self.K_AccV*c_V + self.K_AccOH*c_OH))*(c_V*self.K_AccV + c_OH*self.K_AccOH))
        return (eqOne, eqTwo)



class GridParameters(object):
    def __init__(self, GrainBoundary, approximateStepLength):
        self.layerArea = GrainBoundary.layerArea
        self.numberOfLayers = GrainBoundary.numberOfLayers
        self.layerLength = GrainBoundary.layerLength
        self.approximateStepLength = approximateStepLength
        self.systemWidth = self.getSystemWidth()
        self.numberOfSteps = self.getNumberOfSteps()
        self.stepLength = self.getStepLength()
        self.stepVolume = self.getStepVolume()
        self.initialPhiList = self.getInitialPhiList()

    def getSystemWidth(self):
        return (self.numberOfLayers - 1)*self.layerLength

    def getNumberOfSteps(self):
        #return int(np.ceil((self.systemWidth/self.approximateStepLength)//2)*2 + 1) #old code
        #return (int(np.ceil((self.layerLength/self.approximateStepLength)//2)*2 + 1))*self.numberOfLayers #funker ikke helt
        return (int(np.ceil((self.layerLength/self.approximateStepLength)//2)*2) + 1)*(self.numberOfLayers - 1) + 1 #denne skal funke - sjekk!

    def getStepLength(self):
        return self.systemWidth/self.numberOfSteps

    def getStepVolume(self):
        return self.layerArea*self.getStepLength()

    def getInitialPhiList(self):
        return np.zeros(self.numberOfSteps)


class ChargeDensity(object):
    def __init__(self, grainBoundaryObj, gridParamsObj, accDoping, T, phiList=0):
        self.grainBoundaryObj = grainBoundaryObj
        self.numberOfLayers = grainBoundaryObj.numberOfLayers
        self.c_V_bulk = grainBoundaryObj.bulkConcentrations.c_V_bulk
        self.c_OH_bulk = grainBoundaryObj.bulkConcentrations.c_OH_bulk
        self.c_Acc_bulk = grainBoundaryObj.bulkConcentrations.c_Acc_bulk
        self.c_AccV_bulk = grainBoundaryObj.bulkConcentrations.c_AccV_bulk
        self.c_AccOH_bulk = grainBoundaryObj.bulkConcentrations.c_AccOH_bulk
        self.segregationEnergyLists = grainBoundaryObj.segregationEnergyLists
        self.segregationEnergies = grainBoundaryObj.segregationEnergies
        self.numberOfSteps = gridParamsObj.getNumberOfSteps()
        self.stepVolume = gridParamsObj.getStepVolume()
        self.layerLength = gridParamsObj.layerLength
        self.stepLength = gridParamsObj.stepLength
        self.systemWidth = gridParamsObj.systemWidth
        self.phiList = phiList
        self.accDoping = accDoping
        self.T = T


    def getChargeDensityListJacobi(self):
        SCLconcentrationsJacobi = self.getSCLconcentrationsJacobi()
        c_V_SCL_list, c_OH_SCL_list, c_Acc_SCL_list, c_AccV_SCL_list = SCLconcentrationsJacobi.c_V_SCL_list, SCLconcentrationsJacobi.c_OH_SCL_list, SCLconcentrationsJacobi.c_Acc_SCL_list, SCLconcentrationsJacobi.c_AccV_SCL_list
        chargeDensityListJacobi = np.zeros(self.numberOfSteps)
        for i in range(0, self.numberOfLayers):
            chargeDensity = (e/self.stepVolume)*(2*c_V_SCL_list[i] + c_OH_SCL_list[i] + c_AccV_SCL_list[i] - c_Acc_SCL_list[i])
            #chargeDensity = (self.layerLength/self.stepLength)*e*(2*c_V_SCL_list[i] + c_OH_SCL_list[i] + c_AccV_SCL_list[i] - c_Acc_SCL_list[i])
            j = int(round((i + 0.5)*len(self.phiList)/self.numberOfLayers))
            chargeDensityListJacobi[j] = chargeDensity
        return chargeDensityListJacobi

    def getSCLconcentrationsJacobi(self):
        c_V_SCL_list, c_OH_SCL_list, c_O_SCL_list, c_Acc_SCL_list, c_AccV_SCL_list, c_AccOH_SCL_list = [], [], [], [], [], []
        for layerIndex in range(0, self.numberOfLayers):
            phiIndex = int(round((layerIndex + 0.5)*len(self.phiList)/(self.numberOfLayers)))
            oxygenSitesPerLayer = len(self.segregationEnergyLists.dE_V_list[layerIndex])
            layerConcentrations = self.getLayerConcentrations(layerIndex, phiIndex, oxygenSitesPerLayer)
            c_V_SCL_list.append(layerConcentrations.c_V_layer)
            c_OH_SCL_list.append(layerConcentrations.c_OH_layer)
            c_Acc_SCL_list.append(layerConcentrations.c_Acc_layer)
            c_AccV_SCL_list.append(layerConcentrations.c_AccV_layer)
            c_AccOH_SCL_list.append(layerConcentrations.c_AccOH_layer)
        SCLconcentrationsJacobi = SCLconcentrations(c_V_SCL_list, c_OH_SCL_list, c_Acc_SCL_list, c_AccV_SCL_list, c_AccOH_SCL_list)
        return SCLconcentrationsJacobi

    def getLayerConcentrations(self, layerIndex, phiIndex, oxygenSitesPerLayer):
        c_V_layer, c_OH_layer, c_Acc_layer, c_AccV_layer, c_AccOH_layer = 0, 0, 0, 0, 0
        for siteIndex in range(0, oxygenSitesPerLayer):
            segregationEnergiesAtOneSite = self.getSegregationEnergiesAtOneSite(layerIndex, siteIndex)
            siteConcentrations = self.getSiteConcentrations(segregationEnergiesAtOneSite, self.phiList[phiIndex])
            c_V_layer = c_V_layer + siteConcentrations.c_V_atOneSite/oxygenSitesPerLayer
            c_OH_layer = c_OH_layer + siteConcentrations.c_OH_atOneSite/oxygenSitesPerLayer
            c_Acc_layer = c_Acc_layer + siteConcentrations.c_Acc_atOneSite/oxygenSitesPerLayer
            c_AccV_layer = c_AccV_layer + siteConcentrations.c_AccV_atOneSite/oxygenSitesPerLayer
            c_AccOH_layer = c_AccOH_layer + siteConcentrations.c_AccOH_atOneSite/oxygenSitesPerLayer
        layerConcentrations = LayerConcentrations(c_V_layer, c_OH_layer, c_Acc_layer, c_AccV_layer, c_AccOH_layer)
        return layerConcentrations

    def getSegregationEnergiesAtOneSite(self, layerIndex, siteIndex):
        dE_V_atOneSite = self.segregationEnergyLists.dE_V_list[layerIndex][siteIndex]
        dE_OH_atOneSite = self.segregationEnergyLists.dE_OH_list[layerIndex][siteIndex]
        dE_AccV_atOneSite = self.segregationEnergyLists.dE_AccV_list[layerIndex][siteIndex]
        dE_AccOH_atOneSite = self.segregationEnergyLists.dE_AccOH_list[layerIndex][siteIndex]
        segregationEnergiesAtOneSite = SegregationEnergiesAtOneSite(dE_V_atOneSite, dE_OH_atOneSite, dE_AccV_atOneSite, dE_AccOH_atOneSite)
        return segregationEnergiesAtOneSite

    def getSiteConcentrations(self, segregationEnergiesAtOneSite, phiAtLayer):
        V_ExpTerm = self.getDefectExponentialTerm(2, segregationEnergiesAtOneSite.dE_V_atOneSite, phiAtLayer)
        OH_ExpTerm = self.getDefectExponentialTerm(1, segregationEnergiesAtOneSite.dE_OH_atOneSite, phiAtLayer)
        AccV_ExpTerm = self.getDefectExponentialTerm(1, segregationEnergiesAtOneSite.dE_AccV_atOneSite, phiAtLayer)
        AccOH_ExpTerm = self.getDefectExponentialTerm(0, segregationEnergiesAtOneSite.dE_AccOH_atOneSite, phiAtLayer)
        #The cluster concentrations are calculated first, as they restrict the number of "free" oxygen sites
        acceptorSiteDenominator = self.accDoping + self.c_AccV_bulk*(AccV_ExpTerm - 1) + self.c_AccOH_bulk*(AccOH_ExpTerm - 1)
        c_AccV_atOneSite = self.accDoping*self.c_AccV_bulk*AccV_ExpTerm/acceptorSiteDenominator
        c_AccOH_atOneSite = self.accDoping*self.c_AccOH_bulk*AccOH_ExpTerm/acceptorSiteDenominator
        c_Acc_atOneSite = self.accDoping - c_AccV_atOneSite - c_AccOH_atOneSite
        #Cluster concentrations are then used in the oxygen site denominator
        oxygenSiteDenominator = 3 - c_AccV_atOneSite - c_AccOH_atOneSite + self.c_V_bulk*(V_ExpTerm - 1) + self.c_OH_bulk*(OH_ExpTerm - 1)
        c_V_atOneSite = (3 - c_AccV_atOneSite - c_AccOH_atOneSite)*self.c_V_bulk*V_ExpTerm/oxygenSiteDenominator
        c_OH_atOneSite = (3 - c_AccV_atOneSite - c_AccOH_atOneSite)*self.c_OH_bulk*OH_ExpTerm/oxygenSiteDenominator
        siteConcentrations = SiteConcentrations(c_V_atOneSite, c_OH_atOneSite, c_Acc_atOneSite, c_AccV_atOneSite, c_AccOH_atOneSite)
        return siteConcentrations

    def getDefectExponentialTerm(self, defectCharge, dE_defect, phiAtX):
        return np.exp(-(dE_defect + defectCharge*phiAtX)/(kBev*self.T))

    def getTotalChargeElectronUnitsJacobi(self):
        return (1/e)*self.stepVolume*sum(self.getChargeDensityListJacobi())

    def getCoreConcentrations(self, phiCore, OsitesInCore):
        c_V_core, c_OH_core, c_Acc_core, c_AccV_core, c_AccOH_core = 0, 0, 0, 0, 0
        for siteIndex in range(0, OsitesInCore):
            segregationEnergiesAtOneSite = SegregationEnergiesAtOneSite(self.segregationEnergies.vacancyCoreSegrationEnergy[siteIndex], self.segregationEnergies.protonCoreSegrationEnergy[siteIndex], self.segregationEnergies.acceptorVacancyCoreSegregationEnergy[siteIndex], self.segregationEnergies.acceptorProtonCoreSegregationEnergy[siteIndex])
            siteConcentrations = self.getSiteConcentrations(segregationEnergiesAtOneSite, phiCore)
            c_V_core = c_V_core + siteConcentrations.c_V_atOneSite/OsitesInCore
            c_OH_core = c_OH_core + siteConcentrations.c_OH_atOneSite/OsitesInCore
            c_Acc_core = c_Acc_core + siteConcentrations.c_Acc_atOneSite/OsitesInCore
            c_AccV_core = c_AccV_core + siteConcentrations.c_AccV_atOneSite/OsitesInCore
            c_AccOH_core = c_AccOH_core + siteConcentrations.c_AccOH_atOneSite/OsitesInCore
        coreConcentrations = CoreConcentrations(c_V_core, c_OH_core, c_Acc_core, c_AccV_core, c_AccOH_core)
        return coreConcentrations

    def getCoreCharge(self, phiCore):
        OsitesInCore = len(self.segregationEnergies.vacancyCoreSegrationEnergy)
        OsitesInUnitCell = 3
        numberOfPlanesParallelToGBinUnitCell = 2
        coreUnitCellVolumeRatio = self.grainBoundaryObj.atomicPlaneVolume/self.grainBoundaryObj.unitCellVolume
        self.coreConcentrations = self.getCoreConcentrations(phiCore, OsitesInCore)
        c_V_core, c_OH_core, c_Acc_core, c_AccV_core = self.coreConcentrations.c_V_core, self.coreConcentrations.c_OH_core, self.coreConcentrations.c_Acc_core, self.coreConcentrations.c_AccV_core
        coreCharge = coreUnitCellVolumeRatio*(numberOfPlanesParallelToGBinUnitCell*OsitesInCore/OsitesInUnitCell)*e*(2*c_V_core + c_OH_core + c_AccV_core)
        return coreCharge

    def getSCLchargeBVP(self, SCLarray, phiSCLlist):
        SCLstepLength = SCLarray[1] - SCLarray[0]
        c_V_SCLlist = self.c_V_bulk*np.exp(-(2*phiSCLlist)/(kBev*self.T))
        c_OH_SCLlist = self.c_OH_bulk*np.exp(-(phiSCLlist)/(kBev*self.T))
        c_Acc_SCLlist = np.full(len(phiSCLlist), self.c_Acc_bulk)
        c_AccV_SCLlist = self.c_AccV_bulk*np.exp(-(phiSCLlist)/(kBev*self.T))
        chargeDensityList = (e/self.grainBoundaryObj.unitCellVolume)*(2*c_V_SCLlist + c_OH_SCLlist + c_AccV_SCLlist - c_Acc_SCLlist)

        #totalSCLcharge = sum(chargeDensityList)*SCLstepLength*(self.grainBoundaryObj.atomicPlaneVolume/self.grainBoundaryObj.layerLength)
        totalSCLcharge = sum(chargeDensityList)*SCLstepLength*self.grainBoundaryObj.layerArea
        return 2*totalSCLcharge

    def getSCLchargeAnalytical(self, SCLlength):
        return -2*(self.c_Acc_bulk/self.grainBoundaryObj.unitCellVolume)*e*SCLlength*self.grainBoundaryObj.layerArea



class BulkConcentrations(object):
    def __init__(self, c_V_bulk, c_OH_bulk, c_O_bulk, c_Acc_bulk, c_AccV_bulk, c_AccOH_bulk):
        self.c_V_bulk = c_V_bulk
        self.c_OH_bulk = c_OH_bulk
        self.c_O_bulk = c_O_bulk
        self.c_Acc_bulk = c_Acc_bulk
        self.c_AccV_bulk = c_AccV_bulk
        self.c_AccOH_bulk = c_AccOH_bulk

class SCLconcentrations(object):
    def __init__(self, c_V_SCL_list, c_OH_SCL_list, c_Acc_SCL_list, c_AccV_SCL_list, c_AccOH_SCL_list):
        self.c_V_SCL_list = c_V_SCL_list
        self.c_OH_SCL_list = c_OH_SCL_list
        self.c_Acc_SCL_list = c_Acc_SCL_list
        self.c_AccV_SCL_list = c_AccV_SCL_list
        self.c_AccOH_SCL_list = c_AccOH_SCL_list

class SegregationEnergyLists(object):
    def __init__(self, dE_V_list, dE_OH_list, dE_AccV_list, dE_AccOH_list):
        self.dE_V_list = dE_V_list
        self.dE_OH_list = dE_OH_list
        self.dE_AccV_list = dE_AccV_list
        self.dE_AccOH_list = dE_AccOH_list

class SiteConcentrations(object):
    def __init__(self, c_V_atOneSite, c_OH_atOneSite, c_Acc_atOneSite, c_AccV_atOneSite, c_AccOH_atOneSite):
        self.c_V_atOneSite = c_V_atOneSite
        self.c_OH_atOneSite = c_OH_atOneSite
        self.c_Acc_atOneSite = c_Acc_atOneSite
        self.c_AccV_atOneSite = c_AccV_atOneSite
        self.c_AccOH_atOneSite = c_AccOH_atOneSite

class SegregationEnergiesAtOneSite(object):
    def __init__(self, dE_V_atOneSite, dE_OH_atOneSite, dE_AccV_atOneSite, dE_AccOH_atOneSite):
        self.dE_V_atOneSite = dE_V_atOneSite
        self.dE_OH_atOneSite = dE_OH_atOneSite
        self.dE_AccV_atOneSite = dE_AccV_atOneSite
        self.dE_AccOH_atOneSite = dE_AccOH_atOneSite

class LayerConcentrations(object):
    def __init__(self, c_V_layer, c_OH_layer, c_Acc_layer, c_AccV_layer, c_AccOH_layer):
        self.c_V_layer = c_V_layer
        self.c_OH_layer = c_OH_layer
        self.c_Acc_layer = c_Acc_layer
        self.c_AccV_layer = c_AccV_layer
        self.c_AccOH_layer = c_AccOH_layer

class CoreConcentrations(object):
    def __init__(self, c_V_core, c_OH_core, c_Acc_core, c_AccV_core, c_AccOH_core):
        self.c_V_core = c_V_core
        self.c_OH_core = c_OH_core
        self.c_Acc_core = c_Acc_core
        self.c_AccV_core = c_AccV_core
        self.c_AccOH_core = c_AccOH_core



if __name__ == '__main__':
    import numpy as np
    from physicalConstants import *
    protonSegregationEnergies = [[-0.59, -0.38], [-0.75], [-0.43, -0.62], [-0.33], [-0.36, 0.34]]
    T = 600 #Temperature (in Kelvin), K
    dH = -0.82 #Hydration enthalpy, meV/K
    dS = -0.92 #Hydration entropy, eV
    pH2O = 0.025 #Water vapor pressure, bar
    accDoping = 0.1 #Dopant concentration, site % and formula unit % (they are the same)
    K_AccV = 1
    K_AccOH = 1
    gb = GrainBoundary('210', 4.2349949378839602e-10, True, protonSegregationEnergies, protonSegregationEnergies, 200e-10, T, dH, dS, pH2O, accDoping, K_AccV, K_AccOH)
