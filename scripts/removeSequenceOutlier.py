'''
script que permite eliminar o remover las secuencias de tamano n que no se encuentren dentro
de la media a estudiar con respecto al tipo de secuencia de interes

Se eliminan todas las secuencias con un largo mayor al promedio mas una desviacion estandar y con un
largo menor al promedio menos una desviacion estandar
'''

import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

#exportar a archivo fasta las secuencias agrupadas
def exportFile(dictSequence, nameFileSequence):

    fileExport = open(nameFileSequence, 'w')

    for element in dictSequence:

        header = ">"+element['id']+" "+ element['description']

        sequence = ""
        for i in range(len(element['sequence'])):
            sequence +=element['sequence'][i]

        fileExport.write(header+"\n")
        fileExport.write(sequence+"\n")

    fileExport.close()

#variables input
inputData = sys.argv[1]
lenData = pd.read_csv(sys.argv[2])
pathOutput = sys.argv[3]

#obtenemos el promedio y la desviacion estandar
meanLen = np.mean(lenData['lenData'])
stdLen = np.std(lenData['lenData'])

#definicion de arreglos para filtro
sequences_down = []
sequences_up = []
sequences_inrange = []

#leemos la data y comenzamos a filtrar los elementos
for record in SeqIO.parse(inputData, "fasta"):

    dictSequence = {'sequence': record.seq, 'id': record.id, 'description': record.description}
    lenSequence = len(record.seq)

    if lenSequence>(meanLen+stdLen):
        sequences_up.append(dictSequence)
    elif lenSequence< (meanLen-stdLen) or lenSequence<=50:#minimo de residuos a considerar en una secuencia
        sequences_down.append(dictSequence)
    else:
        sequences_inrange.append(dictSequence)

#exportamos los archivos
exportFile(sequences_up, pathOutput+"removeOver.fasta")
exportFile(sequences_down, pathOutput+"removeUnder.fasta")
exportFile(sequences_inrange, pathOutput+"sequencesInRange.fasta")
