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
def searchPatternClass(sequence, patternClass):


    responseClass = patternClass.findall(str(sequence))

    if len(responseClass)>0:
        return 1
    else:
        return 0

#recibimos los parametros
dataSequence = sys.argv[1]
pathResponse = sys.argv[2]

elementsClassI = []#secuencias clase I
elementsClassII = []#secuencias clase II
otherData = []#no se puede clasificar en algun patron de los reportados
bothData = []#presenta los dos patrones...

#formo el patron de la expresion regular para la clase I
patternClassI = re.compile('C[A-Z]{5,7}CC[A-Z]{19,39}C[A-Z]{8,23}C[A-Z]{5}CC[A-Z]{6,18}C')

#formo el patron de la expresion regular para la clase II
patternClassII = re.compile('C[A-Z]{9,10}CC[A-Z]{11}C[A-Z]{16}C[A-Z]{8,9}CC[A-Z]{10}C')

#hacemos la lectura de las secuencias
for record in SeqIO.parse(dataSequence, "fasta"):

    dictSequence = {'sequence': record.seq, 'id': record.id, 'description': record.description}

    #evaluamos si es clase 1
    response1 = searchPatternClass(record.seq, patternClassI)
    response2 = searchPatternClass(record.seq, patternClassII)

    if response1 == 1:#class I
        elementsClassI.append(dictSequence)
    if response2 == 1:#class II
        elementsClassII.append(dictSequence)

    if response1 == 1 and response2 == 1:#ambas
        bothData.append(dictSequence)

    if response1 == 0 and response2 == 0:#ninguna
        otherData.append(dictSequence)

#exportamos los files
exportFile(elementsClassI, pathResponse+"hydrophobinsClassI.fasta")
exportFile(elementsClassII, pathResponse+"hydrophobinsClassII.fasta")
exportFile(otherData, pathResponse+"otherHydrophobins.fasta")
exportFile(bothData, pathResponse+"bothDataClass.fasta")
