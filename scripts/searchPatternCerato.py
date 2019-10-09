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

def searchPatternCerato(sequence):

    #Analizamos una parte del patron
    patternCerato1 = re.compile('CS[A-Z]{1,3}G[A-Z]{1,3}G')
    patternCerato2 = re.compile('CG[A-Z]{1,3}C')

    response1 = patternCerato1.findall(str(sequence))
    response2 = patternCerato2.findall(str(sequence))

    if len(response1) >0 and len(response2):
        return 1
    else:
        return 0

#recibimos los parametros
dataSequence = sys.argv[1]
pathResponse = sys.argv[2]

elementWithPattern = []

#hacemos la lectura de las secuencias
for record in SeqIO.parse(dataSequence, "fasta"):
    responseClass = searchPatternCerato(record.seq)

    if responseClass == 1:
        dictSequence = {'sequence': record.seq, 'id': record.id, 'description': record.description}
        elementWithPattern.append(dictSequence)


#exportamos los files
exportFile(elementWithPattern, pathResponse+"CeratoProteinsCandidate.fasta")
