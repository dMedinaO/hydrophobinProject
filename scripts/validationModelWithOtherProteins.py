'''
script que recibe un conjunto de datos de espectros de una proteina de interes X y diferentes espectros de
proteinas Y, en este caso: el resto de las generadas en este analisis.

Aplica el metodo de entrenamiento extrayendo un conjunto de testeo de la proteina de interes y formando un conjunto
de testeo de las proteinas restantes para armar clases: es proteina, no es proteina y testear el modelo con estos
elementos con el fin de evaluar si permite solo distinguir proteinas del tipo a analizar.
'''

import sys
from scipy import stats
import pandas as pd
import os
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from sklearn.utils.multiclass import unique_labels
from sklearn.metrics import classification_report
import random

def plot_confusion_matrix(y_true, y_pred, nameFig, normalize=False,title=None,cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = ['Is Protein', 'Is Not Protein']
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()

    plt.savefig(nameFig)
    return ax

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

#funcion que permite extraer espectros aleatorios desde la matriz
def getRandomSpectrum(matrixData, numberElements):

    matrixTesting = []

    for i in range(numberElements):
        index= random.randrange(len(matrixData))
        matrixTesting.append(matrixData[index])

    return matrixTesting

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

#funcion que permite evaluar si un espectro se encuentra en un intervalo de confianza con respecto a las zonas planteadas
def checkClfNewSpectrum(matrixTesting, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax):

    arrayCLF = []

    for element in matrixTesting:
        rowZone = []
        for i in range(len(element)):
            rowZone.append(defineZoneByPoint(element[i], espectralMin[i], espectralMax[i], confidenceIntervalMin[i], confidenceIntervalMax[i]))

        countZ1, countZ2, countZ3, countZ4, countZ5 = (countElementZone(rowZone))

        #clasificamos de la clase, si y solo si el contZ3 es igual al largo del espectro, esto indica que todos los elementos se encuentran dentro
        #del intervalo de confianza del espectro
        if countZ3 >=95:
            arrayCLF.append(0)
        else:
            arrayCLF.append(1)

    return arrayCLF

#obtenemos los nombres de archivos y el path desde la linea de comando
spectrumModel = sys.argv[1]
spectrumTest1 = sys.argv[2]
#spectrumTest2 = sys.argv[3]
#spectrumTest3 = sys.argv[4]
pathOutout = sys.argv[3]
numberTesting = int(sys.argv[4])#numero de secuencias a extraer desde la matriz

#obtenemos las matrices
matrixModel = readDocs(spectrumModel)
matrixTest1 = readDocs(spectrumTest1)
#matrixTest2 = readDocs(spectrumTest2)
#matrixTest3 = readDocs(spectrumTest3)

#obtenemos los valores del intervalo de confianza
espectralMin, espectralMax, espectralMean, confidenceIntervalMin, confidenceIntervalMax = getPropertiesForSpectrum(matrixModel)

testingModel = getRandomSpectrum(matrixModel, numberTesting)
testing1 = getRandomSpectrum(matrixTest1, numberTesting)
#testing2 = getRandomSpectrum(matrixTest2, numberTesting)
#testing3 = getRandomSpectrum(matrixTest3, numberTesting)

lenTotal = len(testingModel)+ len(testing1)

#formamos el arreglo de clases "reales"
realClass = []
for i in range(lenTotal):
    if i <numberTesting:
        realClass.append(0)#si es la proteina
    else:
        realClass.append(1)#no es la proteina

#comenzamos a testear los modelos
predictModel = checkClfNewSpectrum(testingModel, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax)
predictTesting1 = checkClfNewSpectrum(testing1, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax)
#predictTesting2 = checkClfNewSpectrum(testing2, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax)
#predictTesting3 = checkClfNewSpectrum(testing3, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax)

#unimos todos los elementos
predicClass = predictModel+predictTesting1

#obtenemos las medidas de desempeno
print "chec performance"
#obtenemos las medidas de desempeno
nameFig = pathOutout+"confusion_matrixTestingModel.png"
plot_confusion_matrix(realClass, predicClass, nameFig, normalize=True, title='Confusion matrix, without normalization')

print classification_report(realClass, predicClass, target_names=['Is Protein', 'Is Not Protein'])
