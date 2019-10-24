'''
script que recibe un conjunto de secuencias reportadas como hidrofobinas y permite clasificarlas
segun la proteina de interes a evaluar, para ello, se genera un vector de clasificacion, junto
con el porcentaje de las zonas asignadas y la clase en base a la mayoritaria
'''
import sys
from scipy import stats
import pandas as pd
import os
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt

#funcion que permite asignar una zona a un punto
def defineZoneByPoint(pointSpectrum, minPoint, maxPoint, icMinPoint, icMaxPoint):

    zone = ''

    if pointSpectrum> maxPoint:
        zone='Z1'
    elif pointSpectrum<= maxPoint and pointSpectrum > icMaxPoint:
        zone='Z2'
    elif pointSpectrum <= icMaxPoint and pointSpectrum>= icMinPoint:
        zone='Z3'
    elif pointSpectrum <icMinPoint and pointSpectrum >= minPoint:
        zone='Z4'
    else:
        zone='Z5'

    return zone

#funcion que permite contar la cantidad de tipos de zona y retorna un arreglo con el contador
def countElementZone(zoneArray):

    countZ1 = 0
    countZ2 = 0
    countZ3 = 0
    countZ4 = 0
    countZ5 = 0

    for zone in zoneArray:
        if zone == 'Z1':
            countZ1+=1
        if zone == 'Z2':
            countZ2+=1
        if zone == 'Z3':
            countZ3+=1
        if zone == 'Z4':
            countZ4+=1
        if zone == 'Z5':
            countZ5+=1

    countZ1 = float(countZ1)*100/len(zoneArray)
    countZ2 = float(countZ2)*100/len(zoneArray)
    countZ3 = float(countZ3)*100/len(zoneArray)
    countZ4 = float(countZ4)*100/len(zoneArray)
    countZ5 = float(countZ5)*100/len(zoneArray)

    return countZ1, countZ2, countZ3, countZ4, countZ5

#funcion que permite castear la data
def castElement(rowData):

    rowDataCast = []

    for element in rowData:
        rowDataCast.append(float(element))
    return rowDataCast

#funcion que permite leer los documentos y obtener las matrices correspondientes
def readDocs(nameDoc):

    matrixSpectrum = []
    fileOpen = open(nameDoc, 'r')
    line = fileOpen.readline()

    while line:
        matrixSpectrum.append(castElement(line.replace("\n", "").split(",")))
        line = fileOpen.readline()
    fileOpen.close()

    return matrixSpectrum

#funcion que permite caracterizar el conjunto de espectros de frecuencia
def getPropertiesForSpectrum(spectrumData):

    espectralMin = []
    espectralMax = []
    espectralMean = []
    confidenceIntervalMax = []
    confidenceIntervalMin = []

    for i in range(len(spectrumData[0])):
        column = []
        for j in range(len(spectrumData)):
            column.append(float(spectrumData[j][i]))

        espectralMin.append(min(column))
        espectralMax.append(max(column))
        espectralMean.append(np.mean(column))

        meandData = np.mean(column)
        stdData = np.std(column)

        conf_int_a = stats.norm.interval(0.95, loc=meandData, scale=stdData)
        confidenceIntervalMax.append(conf_int_a[1])
        confidenceIntervalMin.append(conf_int_a[0])

    #retornamos los valores correspondientes
    return espectralMin, espectralMax, espectralMean, confidenceIntervalMin, confidenceIntervalMax

#funcion que permite evaluar si un espectro se encuentra en un intervalo de confianza con respecto a las zonas planteadas
def checkClfNewSpectrum(matrixClf, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax):

    arrayClass = []
    arrayZone1 = []
    arrayZone2 = []
    arrayZone3 = []
    arrayZone4 = []
    arrayZone5 = []

    for element in matrixClf:
        rowZone = []
        for i in range(len(element)):
            rowZone.append(defineZoneByPoint(element[i], espectralMin[i], espectralMax[i], confidenceIntervalMin[i], confidenceIntervalMax[i]))

        countZ1, countZ2, countZ3, countZ4, countZ5 = (countElementZone(rowZone))

        #clasificamos de la clase, si y solo si el contZ3 es igual al largo del espectro, esto indica que todos los elementos se encuentran dentro
        #del intervalo de confianza del espectro
        if countZ3 >=100:
            arrayClass.append(0)
        else:
            arrayClass.append(1)
        arrayZone1.append(countZ1)
        arrayZone2.append(countZ2)
        arrayZone3.append(countZ3)
        arrayZone4.append(countZ4)
        arrayZone5.append(countZ5)

    return arrayClass, arrayZone1, arrayZone2, arrayZone3, arrayZone4, arrayZone5

#recibimos la data
spectralSequenceModel = sys.argv[1]#archivo de espectros
inputSpectral = sys.argv[2]#espectros a clasificar
idSequenceData = pd.read_csv(sys.argv[3])#id de secuencias para asociar

spectrumInput = readDocs(spectralSequenceModel)
spectrumClf = readDocs(inputSpectral)

espectralMin, espectralMax, espectralMean, confidenceIntervalMin, confidenceIntervalMax = getPropertiesForSpectrum(spectrumInput)

#obtenemos las clasificaciones
arrayClass, arrayZone1, arrayZone2, arrayZone3, arrayZone4, arrayZone5 = checkClfNewSpectrum(spectrumClf, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax)

#agregamos los elementos al dataFrame
idSequenceData['class'] = arrayClass[:len(arrayClass)-1]
idSequenceData['zone1'] = arrayZone1[:len(arrayClass)-1]
idSequenceData['zone2'] = arrayZone2[:len(arrayClass)-1]
idSequenceData['zone3'] = arrayZone3[:len(arrayClass)-1]
idSequenceData['zone4'] = arrayZone4[:len(arrayClass)-1]
idSequenceData['zone5'] = arrayZone5[:len(arrayClass)-1]

#exportamos el dataFrame
idSequenceData.to_csv(sys.argv[3], index=False)
