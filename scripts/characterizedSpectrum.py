'''
script que recibe una matriz con todos los espectros de frecuencia obtenidos a partir de la digitalizacion
de las propiedades aplicadas en la codificacion de las secuencias y permite caracterizar los espectros, es decir,
obtener propiedades estadisticas del espectro, formar un espectro minimo, maximo y propmedio para luego graficarlo,
genera un espectro multiplicando todos los elementos y sera el espectro caracteristico de la data
'''

import sys
import numpy as np
from scipy import stats
import pandas as pd

def createConsensusPoint(columnData):

    pointData = 1.0

    for element in columnData:
        pointData = float(pointData)*float(element)

    return pointData

espectralDoc = sys.argv[1]#archivo de digitalizacion
pathOutput = sys.argv[2]

#formamos la matriz con la data obtenida desde los espectros
espectralMatrix = []

espectralFile = open(espectralDoc, 'r')
line = espectralFile.readline()

while line:
    row = line.replace("\n", "").split(",")
    espectralMatrix.append(row)
    line = espectralFile.readline()
espectralFile.close()

#obtenemos las caracteristicas del espectro
espectralMin = []
espectralMax = []
espectralMean = []
espectralConsensus = []
confidenceIntervalMax = []
confidenceIntervalMin = []

for i in range(len(espectralMatrix[0])):
    column = []
    for j in range(len(espectralMatrix)):
        column.append(float(espectralMatrix[j][i]))

    espectralMin.append(min(column))
    espectralMax.append(max(column))
    espectralMean.append(np.mean(column))
    espectralConsensus.append(createConsensusPoint(column))

    meandData = np.mean(column)
    stdData = np.std(column)

    conf_int_a = stats.norm.interval(0.95, loc=meandData, scale=stdData)
    confidenceIntervalMax.append(conf_int_a[1])
    confidenceIntervalMin.append(conf_int_a[0])

#generamos una matriz y exportamos la data
matrixValues = [espectralMin, espectralMax, espectralMean, confidenceIntervalMin, confidenceIntervalMax]

dataFrame = pd.DataFrame(matrixValues)
dataFrame.to_csv(pathOutput+"spectrums_properties.csv", index=False)
