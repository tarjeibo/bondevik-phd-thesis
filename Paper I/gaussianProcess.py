from importedLibraries import *

class PredictData(object):
    def __init__(self, importedData, gprParams, energyFilter, indexesToCalculate):
        self.DFTdataDict = importedData.DFTdataDict
        self.descriptorMatrix = importedData.descriptorMatrix
        self.importedData = importedData
        self.gprParams = gprParams
        self.energyFilter = energyFilter
        self.indexesToCalculate = indexesToCalculate
        self.trainGP()
        self.predictGP()

    def trainGP(self):
        self.fetchDataToTrainOn()
        gp_kernel = self.gprParams["gp_kernel"]
        self.gpr = GaussianProcessRegressor(kernel=gp_kernel, n_restarts_optimizer=self.gprParams["n_restarts_optimizer"], alpha=self.gprParams["gpr_noise_param"])
        self.gpr.fit(self.X, self.y)

    def predictGP(self):
        self.y_gpr, self.y_std = self.gpr.predict(self.tempDescriptorMatrix, return_std=True)

    def fetchDataToTrainOn(self):
        self.X = self.descriptorMatrix[self.indexesToCalculate, :]
        self.y = np.take(self.DFTdataDict["energy"], self.indexesToCalculate)
        self.filterAllData()
        self.scaleAllData()

    def filterAllData(self):
        self.trainingIndexesToFilter = []
        self.trainingIndexesToNotFilter = []
        for idx in range(0, len(self.y)):
            if self.y[idx] > min(self.y) + self.energyFilter:
                self.trainingIndexesToFilter.append(idx)
            else:
                self.trainingIndexesToNotFilter.append(idx)
        self.y[self.trainingIndexesToFilter] = min(self.y) + self.energyFilter
        self.tempDescriptorMatrix = np.copy(self.descriptorMatrix)

        # if len(self.trainingIndexesToFilter) > 0:
        #     self.tempDescriptorMatrix[:, 0] = self.descriptorMatrix[:, 0]
        #     for columnIdx in range(1, len(self.X[0])):
        #         maxDescriptor = np.amax(np.take(self.X[:, columnIdx], self.trainingIndexesToNotFilter))
        #         for rowIdx in self.trainingIndexesToFilter:
        #             self.X[rowIdx, columnIdx] = maxDescriptor #Here, we filter the self.X we are about to train
        #         tempDescriptorColumn = []
        #         for rowIdx in range(0, len(self.descriptorMatrix[:, columnIdx])):
        #             if self.descriptorMatrix[rowIdx, columnIdx] > maxDescriptor:
        #                 tempDescriptorColumn.append(maxDescriptor)
        #             else:
        #                 tempDescriptorColumn.append(self.descriptorMatrix[rowIdx, columnIdx])
        #         self.tempDescriptorMatrix[:, columnIdx] = np.asarray(tempDescriptorColumn)

    def scaleAllData(self):
        self.trainingDataScalerDict = {}
        for descriptorIdx in range(0, len(self.X[0])):
            trainingDescriptorColumn = self.X[:, descriptorIdx].reshape(-1, 1)
            descriptorScaler = StandardScaler()
            descriptorScaler.fit(trainingDescriptorColumn)
            self.X[:, descriptorIdx] = descriptorScaler.transform(trainingDescriptorColumn)[:,0]

            completeDescriptorColumn = self.tempDescriptorMatrix[:, descriptorIdx].reshape(-1, 1)
            self.tempDescriptorMatrix[:, descriptorIdx] = descriptorScaler.transform(completeDescriptorColumn)[:,0]

            self.trainingDataScalerDict[descriptorIdx] = descriptorScaler
        self.y = self.y.reshape(-1, 1)
        self.energyScaler = StandardScaler()
        self.energyScaler.fit(self.y)
        self.y = self.energyScaler.transform(self.y)
        self.y = self.y[:,0]

if __name__ == '__main__':
    # print(2/3)
    # a = [[2, 3, 3], [2,3, 5]]
    # a = np.asarray(a)
    # print(len(a))
    # print(len(a[0]))
    # print(np.zeros((2,4)))
    print(np.random.choice(640, 10, replace=False).tolist())
