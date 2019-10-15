'''
script que recibe un archivo csv con el largo de las secuencias analizadas, ademas del archivo multifasta de secuencias
a analizar y como output genera un archivo csv con el ID de la secuencia y si es perteneciente al umbral o no.
# NOTE: Los ID de los elementos trabajan con el identificador de la secuencia. La idea es ir agregando columnas al conjunto de
datos con respecto a cada filtro.
'''

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy import stats

def evalSequence(sequence, conf_int_a):

    response=1
    if len(sequence)>=conf_int_a[0] and len(sequence)<=conf_int_a[1]:
        response=0
    return response

#recibimos el archivo de secuencias y el archivo de tamanos de secuencia
lenSequence = pd.read_csv(sys.argv[2])
sequenceAnalysis = sys.argv[1]
pathResponse = sys.argv[3]

#obtenemos los estadisticos y generamos el intervalo de confianza
meandData = np.mean(lenSequence['lenData'])
stdData = np.std(lenSequence['lenData'])

conf90 = stats.norm.interval(0.90, loc=meandData, scale=stdData)
conf95 = stats.norm.interval(0.95, loc=meandData, scale=stdData)
conf99 = stats.norm.interval(0.99, loc=meandData, scale=stdData)

print conf90
print conf95
print conf99

matrixResponse = []
cont90 = 0
cont95 = 0
cont99 = 0

#recorremos las secuencias y
for record in SeqIO.parse(sequenceAnalysis, "fasta"):

    response90 = evalSequence(record.seq, conf90)
    response95 = evalSequence(record.seq, conf95)
    response99 = evalSequence(record.seq, conf99)

    if response99==0:
        cont99+=1

    if response95==0:
        cont95+=1

    if response90==0:
        cont90+=1

    #formamos la fila
    row = [record.id, response90, response95, response99]
    matrixResponse.append(row)

#formamos un dataframe y exportamos la data
dataframe = pd.DataFrame(matrixResponse, columns=['ID', 'conf90', 'conf95', 'conf99'])
dataframe.to_csv(pathResponse+"analysisPatternData.csv", index=False)

print cont90
print cont95
print cont99
