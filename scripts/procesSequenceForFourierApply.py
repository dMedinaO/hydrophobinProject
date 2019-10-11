'''
script que recibe un conjunto de secuencias y permite poder obtener 6 matrices que reflejan la codificacion y zero-padding de la data para
una posterior digitalizacion de la muestra y obtener los espectros de frecuencia
'''

import sys
from Bio import SeqIO
import numpy as np
import pandas as pd

def exportMatrixData(matrixData, nameDoc):

    dataFrame = pd.DataFrame(matrixData)
    dataFrame.to_csv(nameDoc, index=False)

def makeZeroPadding(rowData):

    lenData = len(rowData)

    pendientes = 512-lenData

    for i in range (pendientes):
        rowData.append(0)

    return rowData

def encodingPropertiesSequence(sequence):

    residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    property1 = [0.61,0.6,0.06,0.46,1.07,0.0,0.47,0.07,0.61,2.22,1.53,1.15,1.18,2.02,1.95,0.05,0.05,2.65,1.88,1.32]#Hydrophobicity
    property2 = [0.046,0.291,0.134,0.105,0.128,0.18,0.151,0.0,0.23,0.186,0.186,0.219,0.221,0.29,0.131,0.062,0.108,0.409,0.298,0.14]#Polarizability
    property3 = [297.0,238.0,236.0,270.0,178.0,185.0,249.0,290.0,277.0,284.0,337.0,224.0,283.0,284.0,222.0,228.0,253.0,282.0,344.0,293.0]#Melting
    property4 = [0.31,-1.01,-0.6,-0.77,1.54,-0.22,-0.64,0.0,0.13,1.8,1.7,-0.99,1.23,1.79,0.72,-0.04,0.26,2.25,0.96,1.22]#Hydrophobic
    property5 = [31.0,124.0,56.0,54.0,55.0,85.0,83.0,3.0,96.0,111.0,111.0,119.0,105.0,132.0,32.5,32.0,61.0,170.0,136.0,84.0]#Volume
    property6 = [-0.5,3.0,0.2,3.0,-1.0,0.2,3.0,0.0,-0.5,-1.8,-1.8,3.0,-1.3,-2.5,0.0,0.3,-0.4,-3.4,-2.3,-1.5]#Hydrophilicity

    encodingDataP1 = []
    encodingDataP2 = []
    encodingDataP3 = []
    encodingDataP4 = []
    encodingDataP5 = []
    encodingDataP6 = []

    for i in range(len(sequence)):

        posData = -1

        for pos in range(len(residues)):#buscamos la posicion del residuo
            if residues[pos] == sequence[i]:
                posData = pos
                break

        #agregamos la propiedad
        encodingDataP1.append(property1[posData])
        encodingDataP2.append(property2[posData])
        encodingDataP3.append(property3[posData])
        encodingDataP4.append(property4[posData])
        encodingDataP5.append(property5[posData])
        encodingDataP6.append(property6[posData])

    #retornamos la data
    return encodingDataP1, encodingDataP2, encodingDataP3, encodingDataP4, encodingDataP5, encodingDataP6

matrixProperty1= []
matrixProperty2= []
matrixProperty3= []
matrixProperty4= []
matrixProperty5= []
matrixProperty6= []

#arreglo del largo de las secuencias
lenElements = []

#datos de consola
inputFile = sys.argv[1]
pathOutput = sys.argv[2]

for record in SeqIO.parse(inputFile, "fasta"):

    lenElements.append(len(record.seq))

meanSequence = np.mean(lenElements)
stdSequence = np.std(lenElements)

#trabajamos solo con las secuencias que se encuentran en un rango con distribucion normal, esto es, valores entre U-2 desv y U+2 desv
sequenceData = []

cortePositivo = meanSequence+1*stdSequence
corteNegativo = meanSequence-1*stdSequence

print max(lenElements)
print min(lenElements)

lenSelected = []
for record in SeqIO.parse(inputFile, "fasta"):

    lenSeq = len(record.seq)

    if lenSeq>corteNegativo:
        if lenSeq < cortePositivo:
            sequenceData.append(record.seq)
            encodingValues = encodingPropertiesSequence(record.seq)

            matrixProperty1.append(encodingValues[0])
            matrixProperty2.append(encodingValues[1])
            matrixProperty3.append(encodingValues[2])
            matrixProperty4.append(encodingValues[3])
            matrixProperty5.append(encodingValues[4])
            matrixProperty6.append(encodingValues[5])
            lenSelected.append(len(record.seq))

#########################################################
print "Corte Positivo: ", cortePositivo
print "Corte Negativo: ", corteNegativo
print max(lenSelected)

matrixProperty1_zeropadding= []
matrixProperty2_zeropadding= []
matrixProperty3_zeropadding= []
matrixProperty4_zeropadding= []
matrixProperty5_zeropadding= []
matrixProperty6_zeropadding= []

for i in range(len(matrixProperty1)):

    matrixProperty1_zeropadding.append(makeZeroPadding(matrixProperty1[i]))
    matrixProperty2_zeropadding.append(makeZeroPadding(matrixProperty2[i]))
    matrixProperty3_zeropadding.append(makeZeroPadding(matrixProperty3[i]))
    matrixProperty4_zeropadding.append(makeZeroPadding(matrixProperty4[i]))
    matrixProperty5_zeropadding.append(makeZeroPadding(matrixProperty5[i]))
    matrixProperty6_zeropadding.append(makeZeroPadding(matrixProperty6[i]))

#exportamos las matrices
exportMatrixData(matrixProperty1_zeropadding, pathOutput+"Hydrophobicity.csv")
exportMatrixData(matrixProperty2_zeropadding, pathOutput+"Polarizability.csv")
exportMatrixData(matrixProperty3_zeropadding, pathOutput+"Melting.csv")
exportMatrixData(matrixProperty4_zeropadding, pathOutput+"Hydrophobic.csv")
exportMatrixData(matrixProperty5_zeropadding, pathOutput+"Volume.csv")
exportMatrixData(matrixProperty6_zeropadding, pathOutput+"Hydrophilicity.csv")
