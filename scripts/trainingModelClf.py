'''
script que permite recibir el conjunto de espectros de frecuencia generadas y a partir de ella obtiene
sub conjuntos de los cuales forma los espectros de frecuencia consenso y en base a estos, genera la clasificacion
de los elementos.

El objetivo general es tomar la data y poder aplicar un sistema de validacion cruzada a los elementos, obtener medidas
de desempeno y formar un intervalo de confianza ponderado con respecto a los puntos a trabajar
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

#funcion que permite evaluar si un espectro se encuentra en un intervalo de confianza con respecto a las zonas planteadas
def checkClfNewSpectrum(spectrumTest, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax):

    rowZone = []
    for i in range(len(spectrumTest[0])):
        rowZone.append(defineZoneByPoint(spectrumTest[0][i], espectralMin[i], espectralMax[i], confidenceIntervalMin[i], confidenceIntervalMax[i]))

    countZ1, countZ2, countZ3, countZ4, countZ5 = (countElementZone(rowZone))

    #clasificamos de la clase, si y solo si el contZ3 es igual al largo del espectro, esto indica que todos los elementos se encuentran dentro
    #del intervalo de confianza del espectro
    if countZ3 >=95:
        return 0
    else:
        return 1

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

        conf_int_a = stats.norm.interval(0.99, loc=meandData, scale=stdData)
        confidenceIntervalMax.append(conf_int_a[1])
        confidenceIntervalMin.append(conf_int_a[0])

    #retornamos los valores correspondientes
    return espectralMin, espectralMax, espectralMean, confidenceIntervalMin, confidenceIntervalMax

#funcion que permite formar el espectro de puntos para entrenar y el espectro para validar
def createDataSetForValidation(matrixSpectrum, index):

    spectrumTest = []
    spectrumTraining = []

    for i in range(len(matrixSpectrum)):
        if i == index:
            spectrumTest.append(matrixSpectrum[i])
        else:
            spectrumTraining.append(matrixSpectrum[i])

    return spectrumTraining, spectrumTest

#funcion que permite castear la data
def castElement(rowData):

    rowDataCast = []

    for element in rowData:
        rowDataCast.append(float(element))
    return rowDataCast

print "GET ELEMENTS INPUT"
#recibimos los datos de interes
espectralValues = sys.argv[1]
pathOutout = sys.argv[2]

#numberSplit = int(sys.argv[3])

print "Read sequence spectrum"
#hacemos la lectura del archivo de espectros
matrixSpectrum = []
fileOpen = open(espectralValues, 'r')
line = fileOpen.readline()

while line:
    matrixSpectrum.append(castElement(line.replace("\n", "").split(",")))
    line = fileOpen.readline()
fileOpen.close()

#con el archivo de datos leido, procedemos a hacer la validaciones correspondientes, para ello, se efectuara un Leave One Out
#esto quiere decir, trabajar con todo el conjunto de datos dejando uno afuera el cual es el que se testeara

print "Training data"
arrayCLF = []
for i in range(len(matrixSpectrum)):

    print "training with sequence: ", i
    #obtenemos el espectro de testeo y el espectro de la data a "entrenar" el modelo
    spectrumTraining, spectrumTest = createDataSetForValidation(matrixSpectrum, i)

    #formamos el intervalo de confianza y los valores de propiedades que permiten caracterizar al espectro
    espectralMin, espectralMax, espectralMean, confidenceIntervalMin, confidenceIntervalMax = getPropertiesForSpectrum(spectrumTraining)

    clfSequence = checkClfNewSpectrum(spectrumTest, espectralMin, espectralMax, confidenceIntervalMin, confidenceIntervalMax)
    arrayCLF.append(clfSequence)

#comenzamos a obtener las medidas de desempeno
arrayReal = []

for i in range(len(arrayCLF)):
    arrayReal.append(0)#es verdad


print "chec performance"
#obtenemos las medidas de desempeno
nameFig = pathOutout+"confusion_matrix.png"
plot_confusion_matrix(arrayReal, arrayCLF, nameFig, normalize=True, title='Confusion matrix, without normalization')

print classification_report(arrayReal, arrayCLF, target_names=['Is Protein', 'Is Not Protein'])
