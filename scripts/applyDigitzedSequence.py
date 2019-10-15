'''
script que permite recibir un conjunto de secuencias, evalua y filtra por largo, las codifica con respecto a las
propiedades de interes, aplica zero-padding y posteriomente ejecuta la rutina matlab que permite obtener el espectro
de frecuencia. Por cada secuencia se crea una carpeta con el ID de la secuencia (de 1 a n solamente) en donde se almacenaran
los resultados para posterior evaluacion de pertenencia al intervalo de confianza de la proteina tipo
'''

import sys
import os
from Bio import SeqIO
import pandas as pd

#funcion que permite hacer el zero-padding
def makeZeroPadding(rowData, maxValue):

    lenData = len(rowData)
    pendientes = maxValue-lenData

    for i in range (pendientes):
        rowData.append(0)
    return rowData

#funcion que permite hacer la codificacion
def encodingPropertiesSequence(sequence):

    residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    property1 = [0.61,0.6,0.06,0.46,1.07,0.0,0.47,0.07,0.61,2.22,1.53,1.15,1.18,2.02,1.95,0.05,0.05,2.65,1.88,1.32]#Hydrophobicity
    property2 = [0.31,-1.01,-0.6,-0.77,1.54,-0.22,-0.64,0.0,0.13,1.8,1.7,-0.99,1.23,1.79,0.72,-0.04,0.26,2.25,0.96,1.22]#Hydrophobic
    property3 = [-0.5,3.0,0.2,3.0,-1.0,0.2,3.0,0.0,-0.5,-1.8,-1.8,3.0,-1.3,-2.5,0.0,0.3,-0.4,-3.4,-2.3,-1.5]#Hydrophilicity

    encodingDataP1 = []
    encodingDataP2 = []
    encodingDataP3 = []

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

    #retornamos la data
    return encodingDataP1, encodingDataP2, encodingDataP3

#recibimos los datos de entrada
sequenceData = sys.argv[1]
pathOutput = sys.argv[2]
maxLenSequence = int(sys.argv[3])

#definicion de matrices y vectores
matrixProperty1 = []
matrixProperty2 = []
matrixProperty3 = []
indexSequence = []

#definicion de valores de propiedade
index = 1
for record in SeqIO.parse(sequenceData, "fasta"):

    if len(record.seq)<=maxLenSequence:

        #aplicamos la codificacion
        encodingDataP1, encodingDataP2, encodingDataP3 = encodingPropertiesSequence(record.seq)

        #agregamos zero-padding
        encodingDataP1ZP =makeZeroPadding(encodingDataP1, maxLenSequence)
        encodingDataP2ZP =makeZeroPadding(encodingDataP2, maxLenSequence)
        encodingDataP3ZP =makeZeroPadding(encodingDataP3, maxLenSequence)

        #agregamos la data a la matriz y al indexSequence
        matrixProperty1.append(encodingDataP1ZP)
        matrixProperty2.append(encodingDataP2ZP)
        matrixProperty3.append(encodingDataP3ZP)
        indexSequence.append(index)
    index+=1

#generamos los data frame y exportamos la data
dataM1 = pd.DataFrame(matrixProperty1)
dataM2 = pd.DataFrame(matrixProperty2)
dataM3 = pd.DataFrame(matrixProperty3)
dataIndex = pd.DataFrame(indexSequence, columns=['indexSequence'])

#exportamos la data
dataM1.to_csv(pathOutput+"Hydrophobicity.csv", index=False)
dataM2.to_csv(pathOutput+"Hydrophobic.csv", index=False)
dataM3.to_csv(pathOutput+"Hydrophilicity.csv", index=False)
dataIndex.to_csv(pathOutput+"indexSequence.csv", index=False)
