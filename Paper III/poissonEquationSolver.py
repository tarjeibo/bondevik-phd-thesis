from importedLibraries import *
from physicalConstants import *
from coordinationComplexConcentration import *

class PoissonEquationJacobiSolver(object):
    #Ref: https://folk.uio.no/jorgeem/download/CSE_example_capacitor_numerical_PDE.pdf
    def __init__(self, chargeDensityObj, jacobiIterations, epsilonR, gridParamsObj, convergenceCriteria=1):
        self.solverType = "Jacobi"
        self.chargeDensityObj = chargeDensityObj
        self.chargeDensityList = chargeDensityObj.getChargeDensityListJacobi()
        self.phiList = chargeDensityObj.phiList
        self.jacobiIterations = jacobiIterations
        self.convergenceCriteria = convergenceCriteria
        self.epsilonR = epsilonR
        self.stepLength = gridParamsObj.stepLength
        self.poissonsEquationSolverJacobi()

    def poissonsEquationSolverJacobi(self):
        N = len(self.phiList)
        epsilon = epsilon0*self.epsilonR
        iterations = 0
        while iterations < self.jacobiIterations:
            phiList_temp = np.copy(self.phiList)
            phiList_temp[0] = 0 #Boundary condition: phi(-zmax)=0
            phiList_temp[N-1] = 0 #Boundary condition: phi(zmax)=0
            for i in range (1,N-1): # Avoid updating the boundaries
                self.phiList[i] = 0.5*(phiList_temp[i+1] + phiList_temp[i-1] + ((self.stepLength**2)/epsilon)*self.chargeDensityList[i])
            self.phiList[0], self.phiList[N-1] = phiList_temp[0], phiList_temp[N-1]
            iterations += 1


class PoissonEquationBVPsolver(object):
    def __init__(self, grainBoundaryObj, gridParamsObj, chargeDensityObj, epsilonR, phiGuess=0.53):
        self.solverType = "BVP"
        self.grainBoundaryObj = grainBoundaryObj
        self.gridParamsObj = gridParamsObj
        self.chargeDensityObj = chargeDensityObj
        self.epsilonR = epsilonR
        self.phiCore = phiGuess
        self.c_V_bulk = self.grainBoundaryObj.bulkConcentrations.c_V_bulk
        self.c_OH_bulk = self.grainBoundaryObj.bulkConcentrations.c_OH_bulk
        self.c_Acc_bulk = self.grainBoundaryObj.bulkConcentrations.c_Acc_bulk
        self.c_AccV_bulk = self.grainBoundaryObj.bulkConcentrations.c_AccV_bulk
        self.c_AccOH_bulk = self.grainBoundaryObj.bulkConcentrations.c_AccOH_bulk
        self.T = self.grainBoundaryObj.T
        self.setSCLarray()
        self.iteratePoissonUntilConvergence()
        self.setMainBVPresults()

    def setSCLarray(self):
        SCLelements = self.gridParamsObj.numberOfSteps//2 + 1
        SCLlength = self.gridParamsObj.systemWidth/2
        self.SCLarray = np.linspace(0, SCLlength, SCLelements)

    def setBVPequationConstants(self):
        self.A = 2*e*self.c_V_bulk/(epsilon0*self.epsilonR)/self.grainBoundaryObj.unitCellVolume
        self.B = e*self.c_OH_bulk/(epsilon0*self.epsilonR)/self.grainBoundaryObj.unitCellVolume
        self.C = e*self.c_Acc_bulk/(epsilon0*self.epsilonR)/self.grainBoundaryObj.unitCellVolume
        self.D = e*self.c_AccV_bulk/(epsilon0*self.epsilonR)/self.grainBoundaryObj.unitCellVolume

    def bvpSolver(self, x, y):
        return np.vstack((y[1], -self.A*np.exp(-(2*y[0])/(kBev*self.T)) - self.B*np.exp(-(y[0])/(kBev*self.T)) + self.C - self.D*np.exp(-(y[0])/(kBev*self.T))))

    def boundaryConditions(self, yStart, yEnd):
        return np.array([yStart[0] - self.phiCore, yEnd[1]])

    def solvePoissonEquationAsBVP(self):
        #Setting initial guess to y = phiGuess, and y' = -1
        y = np.zeros((2, self.SCLarray.size))
        y[0][0] = self.phiCore
        y[1][0] = -1
        res = solve_bvp(self.bvpSolver, self.boundaryConditions, self.SCLarray, y)
        self.phiSCLlist = res.sol(self.SCLarray)[0]

    def iteratePoissonUntilConvergence(self):
        self.setBVPequationConstants()
        self.SCLcharge, self.coreCharge = 1, 1 #dummy values to enter while loop
        while self.coreCharge + self.SCLcharge > 0:
            self.solvePoissonEquationAsBVP()
            self.coreCharge = self.chargeDensityObj.getCoreCharge(self.phiCore)
            self.SCLcharge = self.chargeDensityObj.getSCLchargeBVP(self.SCLarray, self.phiSCLlist)
            self.phiCore = self.phiCore + 0.01

    def setMainBVPresults(self):
        SCLstepLength = self.SCLarray[1] - self.SCLarray[0]
        phiSCLlist = self.phiSCLlist
        c_V_SCLlist = self.c_V_bulk*np.exp(-(2*phiSCLlist)/(kBev*self.T))
        c_OH_SCLlist = self.c_OH_bulk*np.exp(-(phiSCLlist)/(kBev*self.T))
        c_Acc_SCLlist = np.full(len(phiSCLlist), self.c_Acc_bulk)
        c_AccV_SCLlist = self.c_AccV_bulk*np.exp(-(phiSCLlist)/(kBev*self.T))
        chargeDensityList = (e/self.grainBoundaryObj.unitCellVolume)*(2*c_V_SCLlist + c_OH_SCLlist + c_AccV_SCLlist - c_Acc_SCLlist)
        self.mainBVPresults = MainBVPresults(chargeDensityList, c_V_SCLlist, c_OH_SCLlist, c_Acc_SCLlist, c_AccV_SCLlist, phiSCLlist)


class PoissonEquationAnalyticalsolver(object):
    def __init__(self, grainBoundaryObj, gridParamsObj, chargeDensityObj, epsilonR, phiGuess=0.39):
        self.solverType = "Analytical"
        self.grainBoundaryObj = grainBoundaryObj
        self.chargeDensityObj = chargeDensityObj
        self.gridParamsObj = gridParamsObj
        self.epsilonR = epsilonR
        self.phiCore = phiGuess
        self.c_Acc_bulk = grainBoundaryObj.bulkConcentrations.c_Acc_bulk

        self.iteratePoissonUntilConvergence()
        self.setMainAnalyticalResults()

    def getSCLlength(self):
        return np.sqrt((2*epsilon0*self.epsilonR*self.phiCore)/(e*self.c_Acc_bulk/self.grainBoundaryObj.unitCellVolume))

    def iteratePoissonUntilConvergence(self):
        self.SCLcharge, self.coreCharge = 0, 1 #dummy values to enter while loop
        while abs(self.SCLcharge) < self.coreCharge:
            self.coreCharge = self.chargeDensityObj.getCoreCharge(self.phiCore)
            self.SCLcharge = self.chargeDensityObj.getSCLchargeAnalytical(self.getSCLlength())
            self.phiCore = self.phiCore + 0.0001

    def setMainAnalyticalResults(self):
        # SCLlength = np.sqrt((2*epsilon0*self.epsilonR*self.phiCore)/(e*self.c_Acc_bulk/self.grainBoundaryObj.unitCellVolume))
        SCLlength = self.getSCLlength()
        xArray = np.linspace(0, SCLlength, self.gridParamsObj.numberOfSteps//2+1)
        phiSCLlist = []
        for x in xArray:
            phiSCLlist.append(self.phiCore*(x/SCLlength - 1)**2)
        self.mainAnalyticalResults = MainAnalyticalResults(SCLlength, phiSCLlist)


class MainBVPresults(object):
    def __init__(self, chargeDensityList, c_V_SCLlist, c_OH_SCLlist, c_Acc_SCLlist, c_AccV_SCLlist, phiSCLlist):
        self.chargeDensityList = chargeDensityList
        self.c_V_SCLlist = c_V_SCLlist
        self.c_OH_SCLlist = c_OH_SCLlist
        self.c_Acc_SCLlist = c_Acc_SCLlist
        self.c_AccV_SCLlist = c_AccV_SCLlist
        self.phiSCLlist = phiSCLlist


class MainAnalyticalResults(object):
    def __init__(self, SCLlength, phiSCLlist):
        self.SCLlength = SCLlength
        self.phiSCLlist = phiSCLlist

if __name__ == '__main__':
    print("noe")
