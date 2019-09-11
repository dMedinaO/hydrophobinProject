'''
script que permite obtener el patron de clase I de las hidrofobinas
Recibe un archivo fasta con las secuencias y genera un archivo fasta con
las secuencias identificadas como clase I, ademas aqu
'''

import sys
from Bio import SeqIO
import re

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

#identificar si es clase I
def searchPatternClass(sequence):

    #formo el patron de la expresion regular para la clase II
    patternClassI = re.compile('C[A-Z]{5,7}CC[A-Z]{19,39}C[A-Z]{8,23}C[A-Z]{5}CC[A-Z]{6,18}C')

    #formo el patron de la expresion regular para la clase II
    patternClassII = re.compile('C[A-Z]{9,10}CC[A-Z]{11}C[A-Z]{16}C[A-Z]{8,9}CC[A-Z]{10}C')


    responseClassI = patternClassI.findall(str(sequence))
    responseClassII = patternClassII.findall(str(sequence))

    classResponse = 0

    if len(responseClassI) >0:
        classResponse=1
    elif len(responseClassII)>0:
        classResponse=2
    else:
        classResponse=-1

    return classResponse

#recibimos los parametros
dataSequence = sys.argv[1]
pathResponse = sys.argv[2]

elementsClassI = []#secuencias clase I
elementsClassII = []#secuencias clase II
otherData = []#no se puede clasificar en algun patron de los reportados

#hacemos la lectura de las secuencias
for record in SeqIO.parse(dataSequence, "fasta"):
    responseClass = searchPatternClass(record.seq)
    dictSequence = {'sequence': record.seq, 'id': record.id, 'description': record.description}

    if responseClass == 1:
        elementsClassI.append(dictSequence)

    elif responseClass == 2:
        elementsClassII.append(dictSequence)
    else:
        otherData.append(dictSequence)

#exportamos los files
exportFile(elementsClassI, pathResponse+"hydrophobinsClassI.fasta")
exportFile(elementsClassII, pathResponse+"hydrophobinsClassII.fasta")
exportFile(otherData, pathResponse+"otherHydrophobins.fasta")
