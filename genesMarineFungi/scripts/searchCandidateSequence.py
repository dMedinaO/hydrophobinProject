'''
script que permite hacer la busqueda de candidatos de hydrophobin en el genoma con respecto a los patrones descritos
en la bibliografia

Input: archivo con secuencias multifasta
Output: archivos en formato fasta con candidatos de clase tipo I y clase tipo II
'''

import sys
from Bio import SeqIO

#funcion que permite buscar el patron correspondiente a las clase II, solo un segmento caracteristico
def searchPatternClassII(sequence):

    #primero: encontrar la posicion del doble enlace
    indexTwoBonds = []
    response = 0
    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            indexTwoBonds.append(i)

    for index in indexTwoBonds:#buscamos si 12 posiciones posteriores existe una C
        try:
            if sequence[index+12] == 'C' or sequence[index+11] == 'C':
                response = 1
        except:
            pass
    return response

#funcion que permite generar un archivo con las secuencias existentes en un arreglo
def createExportFile(arraySequence, fileExport):

    fileOpen = open(fileExport, 'w')
    index = 0
    for sequence in arraySequence:

        sequenceData = ""
        for i in range(len(sequence)):
            sequenceData+=sequence[i]

        fileOpen.write(">"+str(index)+"\n")
        fileOpen.write(sequenceData)
        fileOpen.write("\n")
        index+=1
    fileOpen.close()

#funcion que permite identificar la diferencia entre las cisteinas
def searchDifferenceCys(sequence):

    indexCys = []

    for i in range(len(sequence)):
        if sequence[i] == 'C':
            indexCys.append(i)

    responseType = []
    for i in range(len(indexCys)-1):
        diff = indexCys[i+1] - indexCys[i]
        try:
            if diff>=5 and diff <=7:
                responseType.append(1)
            elif diff>=9 and diff<=10:
                responseType.append(2)
            else:
                responseType.append(0)
        except:
            responseType.append(0)
            pass
    #con respecto a los posibles tipos evaluo cual es la opcion: 0-> nada, 1-> class I, 2-> class II, 3-> ambas
    contData = 0
    exist1 = 0
    exist2 = 0

    for element in responseType:
        if element == 0:
            contData+=1
        elif element == 1:
            exist1 = 1
        else:
            exist2 = 1

    if len(responseType) == contData:#de ninguna
        return 0
    elif exist1 == 1 and exist2 == 1:#de ambas
        return 3
    elif exist1 == 1 and exist2 == 0:#solo de clase I
        return 1
    else:#solo de clase II
        return 2

#funcion que permite identificar las secuencias que tienen dobles enlaces de cisteinas
def searchSequenceWithTwoCisteins(sequence):

    cont=0
    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            cont+=1

    if cont>=2:
        return 0#posee doble enlace de cisteinas
    else:
        return 1#no posee doble enlace de cisteinas

sequenceCandidateTwoCys = []

#hacemos la lectura del archivo de secuencias
for record in SeqIO.parse(sys.argv[1], "fasta"):

    sequence = record.seq
    if searchSequenceWithTwoCisteins(sequence) == 0:
        sequenceCandidateTwoCys.append(sequence)

createExportFile(sequenceCandidateTwoCys, "sequencesTwoCysteine.fasta")
