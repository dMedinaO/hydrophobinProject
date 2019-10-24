'''
script que permite poder obtener las secuencias de las hidrofobinas unknown en base a la clase, forma un nuevo fasta
con la data que obtiene
'''

import sys
import pandas as pd
from Bio import SeqIO

#funcion que permite generar un archivo con las secuencias existentes en un arreglo
def createExportFile(arraySequence, idSequences, descriptionSequence, fileExport):

    fileOpen = open(fileExport, 'w')
    for i in range(len(arraySequence)):

        sequenceData = ""
        for j in range(len(arraySequence[i])):
            sequenceData+=arraySequence[i][j]

        header = idSequences[i]+" "+ descriptionSequence[i]
        fileOpen.write(">"+header+"\n")
        fileOpen.write(sequenceData)
        fileOpen.write("\n")
    fileOpen.close()

#recibimos la data
sequenceFile = sys.argv[1]
classFile = pd.read_csv(sys.argv[2])
pathResponse = sys.argv[3]

#hacemos la lectura y obtenemos la secuencia
arraySequence = []
arrayDescription = []
arrayIds = []

#hacemos la lectura del archivo de secuencias db1
for record in SeqIO.parse(sequenceFile, "fasta"):

    for element in classFile['indexSequence']:

        if record.id == element:
            arraySequence.append(record.seq)
            arrayIds.append(record.id)
            arrayDescription.append(record.description)
            break

exportData = pathResponse+"sequenceInclass.fasta"
createExportFile(arraySequence, arrayIds, arrayDescription, exportData)
