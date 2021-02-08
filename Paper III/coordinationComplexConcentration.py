#Module to calculate (Acc-OH)-concentrations
from importedLibraries import *
from physicalConstants import *

class CoordinationComplexConcentration(object):
    def __init__(self, grainBoundaryObj, dH_AccV, dH_AccOH, phi, dE_seg_V, dE_seg_OH):
        self.T = grainBoundaryObj.T
        self.dH_hydr = grainBoundaryObj.dH
        self.dS_hydr = grainBoundaryObj.dS
        self.pH2O = grainBoundaryObj.pH2O
        self.c_Acc = grainBoundaryObj.c_Acc
        self.c_Acc_initial = grainBoundaryObj.c_Acc
        self.dH_AccV = dH_AccV
        self.dH_AccOH = dH_AccOH
        self.phi = phi
        self.dE_seg_V = dE_seg_V
        self.dE_seg_OH = dE_seg_OH

        self.kappa = self.getKappa()
        self.setK_AccV()
        self.setK_AccOH()
        self.setAllPseudoConcentrations()
        self.setConcentrationsConsideringSiteRestriction()
        # self.printAllConcentrations()
        # self.printConcentrationsWithNoCC()

    def setK_AccV(self):
        self.K_AccV = np.exp(-self.dH_AccV/(kBev*self.T))

    def setK_AccOH(self):
        self.K_AccOH = np.exp(-self.dH_AccOH/(kBev*self.T))

    def setAllPseudoConcentrations(self):
        self.c_V_pseudo, self.c_OH_pseudo = fsolve(self.coordinationComplexSolver, (0.5, 0.5))
        self.c_Acc = self.getAcceptorConcentration(self.c_V_pseudo, self.c_OH_pseudo)
        self.c_AccV = self.K_AccV*self.c_Acc*self.c_V_pseudo
        self.c_AccOH = self.K_AccOH*self.c_Acc*self.c_OH_pseudo

    def setConcentrationsConsideringSiteRestriction(self):
        self.c_V = 3*self.c_V_pseudo*np.exp(-(self.dE_seg_V + 2*self.phi)/(kBev*self.T))/(3 + self.c_V_pseudo*(np.exp(-(self.dE_seg_V + 2*self.phi)/(kBev*self.T)) - 1) + self.c_OH_pseudo*(np.exp(-(self.dE_seg_OH + self.phi)/(kBev*self.T)) - 1))
        self.c_OH = 3*self.c_OH_pseudo*np.exp(-(self.dE_seg_OH + self.phi)/(kBev*self.T))/(3 + self.c_V_pseudo*(np.exp(-(self.dE_seg_V + 2*self.phi)/(kBev*self.T)) - 1) + self.c_OH_pseudo*(np.exp(-(self.dE_seg_OH + self.phi)/(kBev*self.T)) - 1))

    def printAllConcentrations(self):
        print("Concentrations with complex formation:")
        print("c_V =", self.c_V)
        print("c_OH =", self.c_OH)
        print("c_Acc =", self.c_Acc)
        print("c_AccV =", self.c_AccV)
        print("c_AccOH =", self.c_AccOH)


    def getAcceptorConcentration(self, c_V, c_OH):
        return self.c_Acc/(1 + self.K_AccV*c_V + self.K_AccOH*c_OH)

    def getKappa(self):
        return np.exp((-self.dH_hydr/(kBev*self.T) + self.dS_hydr*1e-3/kBev))*self.pH2O

    def lindmanEquationsToSolve(self, variables):
        c_V_noCC, c_OH_noCC = variables
        return self.kappa*c_V_noCC*(3 - c_V_noCC - c_OH_noCC) - c_OH_noCC**2, 2*c_V_noCC + c_OH_noCC - self.c_Acc_initial

    def coordinationComplexSolver(self, variables):
        c_V, c_OH = variables
        rightSideEq = self.c_Acc*(1 - self.K_AccV*c_V)/(1 + self.K_AccV*c_V + self.K_AccOH*c_OH)
        return self.kappa*c_V*(3 - c_V - c_OH) - c_OH**2, 2*c_V + c_OH - rightSideEq

    def printConcentrationsWithNoCC(self):
        c_V_noCC, c_OH_noCC = fsolve(self.lindmanEquationsToSolve, (0.5, 0.5))
        print("\n")
        print("Concentrations without complex formation:")
        print("c_V =", c_V_noCC)
        print("c_OH =", c_OH_noCC)

if __name__ == '__main__':
    from scipy.optimize import fsolve, newton, root, broyden2
    from physicalConstants import *
    import math
    import numpy as np
    #
    # T = 1000
    # pH2O = 0.025
    # dH_hydr = -0.82
    # dS_hydr = -0.92
    # c_Acc = 0.1
    # #dH_AccV = 0.02
    # #dH_AccOH = -0.53
    # dH_AccV = 0
    # dH_AccOH = 0
    # phi = 0
    # dE_seg_V = 0
    # dE_seg_OH = 0
    # ccc1 = CoordinationComplexConcentration(T, pH2O, dH_hydr, dS_hydr, c_Acc, dH_AccV, dH_AccOH, phi, dE_seg_V, dE_seg_OH)


    #

    def testFunksjon():
        c_V, c_OH = variables
        rightSideEq = self.c_Acc*(1 - self.K_AccV*c_V)/(1 + self.K_AccV*c_V + self.K_AccOH*c_OH)
        return self.kappa*c_V*(3 - c_V - c_OH) - c_OH**2, 2*c_V + c_OH - rightSideEq


    def equations(variables):
        c_V, c_OH = variables
        kappa = 178
        K_AccV = 1
        K_AccOH = 1
        D = 0.1
        eqOne = 2*c_V + c_OH - D*(1 - c_V*K_AccV)/(1 + K_AccV*c_V + K_AccOH*c_OH)
        eqTwo = c_OH**2 + (c_V*kappa)*(-3 + c_V + c_OH + (D/(1 + K_AccV*c_V + K_AccOH*c_OH))*(c_V*K_AccV + c_OH*K_AccOH))
        return (eqOne, eqTwo)


    c_V, c_OH =  fsolve(equations, (1, 1))
    print(c_V, c_OH)

    #
    #
    # def equations(variables):
    #     c_V, c_OH, c_O, c_Acc, c_AccV, c_AccOH = variables
    #     kappa = 178
    #     K_AccV = 1195
    #     K_AccOH = 13
    #     D = 0.1
    #     eqOne = c_V + c_OH + c_O + c_AccV + c_AccOH - 3
    #     eqTwo = 2*c_V + c_OH + c_AccV - c_Acc
    #     eqThree = kappa - (c_OH**2)/(c_V*c_O)
    #     eqFour = K_AccV - c_AccV/(c_Acc*c_V)
    #     eqFive = K_AccOH - c_AccOH/(c_Acc*c_OH)
    #     eqSix = c_Acc + c_AccV + c_AccOH - D
    #     return (eqOne, eqTwo, eqThree, eqFour, eqFive, eqSix)
    #
    #
    # c_V, c_OH, c_O, c_Acc, c_AccV, c_AccOH =  fsolve(equations, (0.01, 0.1, 0.5, 0.1, 0.01, 0.01))
    # print(c_V, c_OH, c_O, c_Acc, c_AccV, c_AccOH)
    #
    # #print(equations((x, y, z)))
    #








#mer plass
