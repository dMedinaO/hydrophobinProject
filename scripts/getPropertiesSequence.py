'''
script que permite recibir un archivo multifasta y obtener propiedades de las secuencias tales como:

1. Largo maximo
2. Largo minimo
3. Largo promedio
4. Desviacion estandar
5. Frecuencia de residuos promedio
'''

import sys
from Bio import SeqIO
import numpy as np
import pandas as pd

#funcion que permite contar los doble cys
def getDobleCys(sequence):

    cont=0

    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            cont+=1

    return cont

#funcion que obtiene las frecuencias
def getFrequenceResidueForSequence(sequence):

    arrayRes = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    arrayFrequence = []

    for residue in arrayRes:
        cont=0
        for i in range(len(sequence)):
            if residue == sequence[i]:
                cont+=1
        cont = float(cont)*100/float(len(sequence))
        arrayFrequence.append(cont)

    return arrayFrequence

lenSequence = []
matrixFrequence = []
dobleCys = []

inputFile = sys.argv[1]
pathOutput = sys.argv[2]

for record in SeqIO.parse(inputFile, "fasta"):

    lenSequence.append(len(record.seq))
    matrixFrequence.append(getFrequenceResidueForSequence(record.seq))
    dobleCys.append(getDobleCys(record.seq))

mean = np.mean(lenSequence)
max = max(lenSequence)
min = min(lenSequence)
std = np.std(lenSequence)

maxDobleCys = np.max(dobleCys)
minDobleCys = np.min(dobleCys)
meanDobleCys = np.mean(dobleCys)
stdDobleCys = np.std(dobleCys)

header = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
#formamos un dataframe con la matriz y el header
dataframe = pd.DataFrame(matrixFrequence, columns=header)

dataframe.to_csv(pathOutput+"frequenceResidues.csv", index=False)

#obtenemos los promedios de cada frecuencia
meandData = []
stdData = []
maxData = []
minData = []

for element in header:
    meandData.append(np.mean(dataframe[element]))
    stdData.append(np.std(dataframe[element]))
    maxData.append(np.max(dataframe[element]))
    minData.append(np.min(dataframe[element]))

matrixStatistics = [meandData, stdData, maxData, minData]

dataFrameStatistics = pd.DataFrame(matrixStatistics, columns=header)
dataFrameStatistics.to_csv(pathOutput+"frequenceResidues_statistics.csv", index=False)

#generamos el archivo resumen
fileOpen = open(pathOutput+"summaryData.txt", 'w')

fileOpen.write("INFORMACION SEQUENCE\n")
fileOpen.write("Largo Promedio: "+str(mean)+"\n")
fileOpen.write("Largo Maximo: "+str(max)+"\n")
fileOpen.write("Largo Minimo: "+str(min)+"\n")
fileOpen.write("Desviacion Largo: "+str(std)+"\n")

fileOpen.write("INFORMACION DOBLE CYS\n")
fileOpen.write("Doble Cys Promedio: "+str(meanDobleCys)+"\n")
fileOpen.write("Doble Cys Maximo: "+str(maxDobleCys)+"\n")
fileOpen.write("Doble Cys Minimo: "+str(minDobleCys)+"\n")
fileOpen.write("Desviacion Doble Cys: "+str(stdDobleCys)+"\n")

fileOpen.close()
