from importedLibraries import *

class ImportAndPreProcessData(object):
    def __init__(self, dataFileName, descriptors, scaleDescriptors=True):
        self.dataFileName = dataFileName
        self.descriptors = descriptors
        self.scaleDescriptors = scaleDescriptors
        self.importDFTdata()
        self.setDescriptorMatrix()

    def importDFTdata(self):
        pbcIndex, aseIndex, xDirect, yDirect, dAngstrom, energy, prePES, bondLengthRatio, nsw3, nsw6, energyCPU, prePESCPU, nsw3CPU, nsw6CPU = \
        [], [], [], [], [], [], [], [], [], [], [], [], [], []
        with open(self.dataFileName) as f:
            reader = csv.reader(f)
            firstLine = True
            for row in reader:
                if firstLine:
                    firstLine = False
                    continue
                pbcIndex.append(float(row[0].split(";")[0]))
                aseIndex.append(float(row[0].split(";")[1]))
                xDirect.append(float(row[0].split(";")[2]))
                yDirect.append(float(row[0].split(";")[3]))
                dAngstrom.append(float(row[0].split(";")[4]))
                energy.append(float(row[0].split(";")[5]))
                prePES.append(float(row[0].split(";")[6]))
                bondLengthRatio.append(float(row[0].split(";")[7]))
                nsw3.append(float(row[0].split(";")[8]))
                nsw6.append(float(row[0].split(";")[9]))
                energyCPU.append(float(row[0].split(";")[10]))
                prePESCPU.append(float(row[0].split(";")[11]))
                nsw3CPU.append(float(row[0].split(";")[12]))
                nsw6CPU.append(float(row[0].split(";")[13]))
        self.DFTdataDict = {
        "pbcIndex": np.asarray(pbcIndex),
        "aseIndex": np.asarray(aseIndex),
        "xDirect": np.asarray(xDirect),
        "yDirect": np.asarray(yDirect),
        "dAngstrom": np.asarray(dAngstrom),
        "energy": np.asarray(energy),
        "prePES": np.asarray(prePES),
        "bondLengthRatio": np.asarray(bondLengthRatio),
        "nsw3": np.asarray(nsw3),
        "nsw6": np.asarray(nsw6),
        "energyCPU": np.asarray(energyCPU),
        "prePESCPU": np.asarray(prePESCPU),
        "nsw3CPU": np.asarray(nsw3CPU),
        "nsw6CPU": np.asarray(nsw6CPU),
        "xDirectUnScaled": np.asarray(xDirect),
        "unScaledEnergy": np.asarray(energy),
        "unScaledLogTransformedEnergy": np.asarray(energy),
        "unScaledPrePES": np.asarray(prePES),
        "diffPES": np.asarray(energy) - np.asarray(prePES)
        }

    def setDescriptorMatrix(self):
        counter = 0
        for descriptor in self.descriptors:
            descriptorColumn = self.DFTdataDict[descriptor]
            if counter == 0:
                self.descriptorMatrix = descriptorColumn
            else:
                self.descriptorMatrix = np.column_stack((self.descriptorMatrix, descriptorColumn))
            counter = counter + 1
        #In the case of only one descriptor, we need the descriptor matrix on 'matrix' form, and not as a column vector
        if counter == 1:
            self.descriptorMatrix = self.descriptorMatrix.reshape(-1, 1)
