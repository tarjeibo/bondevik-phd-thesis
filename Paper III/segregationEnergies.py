#File with segregation energies

def getSegregationEnergies(hkl):

    if hkl == '210':
        vacancySegregationEnergies = [[0.00, 0.16], [0.01], [0.15, -0.04], [-0.14], [-0.08, 0.08], [-0.16], [0.10, -1.41], [-1.02], [-1.02, 0.63], [-1.41], [-1.83, -0.42], [-0.52], [1.79, -0.33]\
        , [-0.11], [0.05, -0.29], [-0.07], [0.03, 0.30], [0.10], [0.21, 0.09], [0.01], [0.00, 0.16]]
        protonSegregationEnergies = [[0.09, 0.00], [-0.03], [0.12, -0.03], [-0.07], [-0.02, 0.12], [-0.18], [0.14, 0.08], [-0.30], [-1.30, 0.52], [-0.71], [-0.71, -0.77], [0.13], [0.21, -0.12]\
        , [-0.13], [0.10, -0.10], [-0.08], [0.01, 0.16], [-0.01], [0.10, 0.00], [-0.03], [0.09, 0.0]]
        acceptorVacancySegregationEnergies = [[0.00, 0.00], [0.00], [-0.17, -0.17], [-0.17], [-0.15, -0.15], [-0.15], [-0.15, -0.15], [-0.15], [-1.69, -1.69], [-1.69], [1.48, 1.48], [-0.05], \
        [-0.05, -0.05], [-0.06], [-0.06, -0.06], [-0.15], [-0.15, -0.15], [-0.15], [-0.15, -0.15], [0.00], [0.00, 0.00]]
        acceptorProtonSegregationEnergies = [[0.00, 0.00], [0.00], [-0.12, -0.12], [-0.12], [-0.09, -0.09], [-0.09], [-0.11, -0.11], [-0.11], [-0.77, -0.77], [-0.77], [0.72, 0.72], [-0.10], \
        [-0.10, -0.10], [-0.10], [-0.10, -0.10], [-0.11], [-0.11, -0.11], [-0.10], [-0.10, -0.10], [0.00], [0.00, 0.00]]
        vacancyCoreSegrationEnergy = [-1.83, -0.42]
        protonCoreSegrationEnergy = [-0.71, -0.77]
        acceptorVacancyCoreSegregationEnergy = [-1.69, -1.69]
        acceptorProtonCoreSegregationEnergy = [-0.77, -0.77]

    elif hkl == '111': #One sided input
        # vacancySegregationEnergies = [[-0.56, -0.56, -0.56], [0, 0, 0], [0, 0, 0], [0.00, 0.00, 0.00]]
        # protonSegregationEnergies = [[-0.74, -0.74, -0.74], [0, 0, 0], [0, 0, 0], [0.00, 0.00, 0.00]]
        vacancySegregationEnergies = [[-0.56, -0.56, -0.56], [0.40, 0.40, 0.40], [0.06, 0.06, 0.06], [0.00, 0.00, 0.00]]
        protonSegregationEnergies = [[-0.74, -0.74, -0.74], [0.17, 0.17, 0.17], [0.10, 0.10, 0.10], [0.00, 0.00, 0.00]]
        acceptorVacancySegregationEnergies = [[0.00, 0.00, 0.00], [-0.08, -0.08, -0.08], [-0.07, -0.07, -0.07], [0.0, 0.0, 0.0]]
        acceptorProtonSegregationEnergies = [[-0.15, -0.15, -0.15], [-0.13, -0.13, -0.13], [-0.10, -0.10, -0.10], [0, 0, 0]]
        vacancyCoreSegrationEnergy = [-0.56, -0.56, -0.56]
        protonCoreSegrationEnergy = [-0.74, -0.74, -0.74]
        acceptorVacancyCoreSegregationEnergy = [0.00, 0.00, 0.00]
        acceptorProtonCoreSegregationEnergy = [-0.15, -0.15, -0.15]
        segregationEnergies = SegregationEnergies(vacancySegregationEnergies, protonSegregationEnergies, acceptorVacancySegregationEnergies, acceptorProtonSegregationEnergies, vacancyCoreSegrationEnergy, protonCoreSegrationEnergy, acceptorVacancyCoreSegregationEnergy, acceptorProtonCoreSegregationEnergy)

    segregationEnergies = SegregationEnergies(vacancySegregationEnergies, protonSegregationEnergies, acceptorVacancySegregationEnergies, acceptorProtonSegregationEnergies, vacancyCoreSegrationEnergy, protonCoreSegrationEnergy, acceptorVacancyCoreSegregationEnergy, acceptorProtonCoreSegregationEnergy)
    return segregationEnergies


class SegregationEnergies(object):
    def __init__(self, vacancySegregationEnergies, protonSegregationEnergies, acceptorVacancySegregationEnergies, acceptorProtonSegregationEnergies, vacancyCoreSegrationEnergy, protonCoreSegrationEnergy, acceptorVacancyCoreSegregationEnergy, acceptorProtonCoreSegregationEnergy):
        self.vacancySegregationEnergies = vacancySegregationEnergies
        self.protonSegregationEnergies = protonSegregationEnergies
        self.acceptorVacancySegregationEnergies = acceptorVacancySegregationEnergies
        self.acceptorProtonSegregationEnergies = acceptorProtonSegregationEnergies
        self.vacancyCoreSegrationEnergy = vacancyCoreSegrationEnergy
        self.protonCoreSegrationEnergy = protonCoreSegrationEnergy
        self.acceptorVacancyCoreSegregationEnergy = acceptorVacancyCoreSegregationEnergy
        self.acceptorProtonCoreSegregationEnergy = acceptorProtonCoreSegregationEnergy





'''Helgee's energies, 210GB:'''
# vacancySegregationEnergies = [[0.39, -0.06], [-0.45], [-0.29, 0.05], [0.00], [-0.47, 0.91]]
# protonSegregationEnergies = [[-0.59, -0.38], [-0.75], [-0.43, -0.62], [-0.33], [-0.36, 0.34]]
#
# vacancyCoreSegrationEnergy = -0.75
# protonCoreSegrationEnergy = -0.45



'''Test case to check Jacobi performance on single plane GB'''
# vacancySegregationEnergies = [[-0.62, -0.62, -0.62], [0.00, 0.00, 0.00], [0.00, 0.00, 0.00]] #Test case, to check how Jacobi performs on single plane GB
# protonSegregationEnergies = [[-0.84, -0.79, -0.75], [0.00, 0.00, 0.00], [0.00, 0.00, 0.00]] #Test case, to check how Jacobi performs on single plane GB
