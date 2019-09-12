'''
script que permite obtener la descripcion de la secuencia, desde la base de datos con informacion
redundante, exporta el mismo archivo pero con la descripcion agregada
'''
import sys
from Bio import SeqIO

def exportFile(sequenceData, descriptionData, idData, fileData):

    fileExport = open(fileData, 'w')

    for i in range(len(sequenceData)):

        header = ">"+idData[i]+" "+descriptionData[i]

        sequence = ""
        for j in range(len(sequenceData[i])):
            sequence+=sequenceData[i][j]

        fileExport.write(header+"\n")
        fileExport.write(sequence+"\n")

    fileExport.close()

classFile = sys.argv[1]
dataBaseFile = sys.argv[2]

sequenceData = []
descriptionData = []
idData = []

for record in SeqIO.parse(classFile, "fasta"):

    description = ""
    print "process sequence ", record.id
    for record2 in SeqIO.parse(dataBaseFile, "fasta"):

        if record2.seq == record.seq:
            description = description+" "+record2.description

    descriptionData.append(description)
    sequenceData.append(record.seq)
    idData.append(record.id)

#exportamos el archivo
exportFile(sequenceData, descriptionData, idData, classFile)
