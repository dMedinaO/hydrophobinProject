'''
script que permite aplicar el filtro de motivos conservados en las secuencias para proteinas del tipo
hidrofobina clase I

Los filtros aplicados son:

1. Numero de cys
2. Numero de cys-cys
3. Patron caracteristico clase I
'''

import sys
from Bio import SeqIO
import pandas as pd
import re

#funcion que permite contar el numero de cys, debe ser mayor o igual a 8
def countNumberCys(sequence):

    contCys=0

    for i in range(len(sequence)):
        if sequence[i] == 'C':
            contCys+=1
    if contCys>=8:#cumple con el minimo
        return 0
    else:#no lo cumple
        return 1

#funcion que permite contar los cys-cys, debe ser mayor o igual a dos
def countNumberCysCys(sequence):

    contCysCys = 0

    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            contCysCys+=1
    if contCysCys>=2:#cumple con el minimo
        return 0
    else:#no lo cumple
        return 1

#funcion que permite evaluar el patron de clase I reportado en la literatura
def evaluatedPattern(sequence):

    patternClassI = re.compile('C[A-Z]{5,7}CC[A-Z]{19,39}C[A-Z]{8,23}C[A-Z]{5}CC[A-Z]{6,18}C')
    responseClass = patternClassI.findall(str(sequence))

    if len(responseClass)>0:#cumple con el patron
        return 0
    else:
        return 1

#recibimos la data de interes
sequenceData = sys.argv[1]
patterAnalysis = pd.read_csv(sys.argv[2])

filterNumberCys = []
cont_filterNumberCys = 0

filterNumberCysCys = []
cont_filterNumberCysCys = 0

filterPattern = []
contPattern = 0

#hacemos la lectura
for record in SeqIO.parse(sequenceData, "fasta"):

    #analizamos el cont cys
    response = countNumberCys(record.seq)
    if response == 0:
        cont_filterNumberCys+=1
    filterNumberCys.append(response)

    #analizamos el cont cys-cys
    response2 = countNumberCysCys(record.seq)
    if response2 == 0:
        cont_filterNumberCysCys+=1
    filterNumberCysCys.append(response2)

    #analizamos el patron reportado
    response3 = evaluatedPattern(record.seq)
    if response3 == 0:
        contPattern+=1
    filterPattern.append(response3)

print cont_filterNumberCys
print cont_filterNumberCysCys
print contPattern

#agregamos la data al dataFrame y exportamos los analisis
patterAnalysis['numberCys'] = filterNumberCys
patterAnalysis['numberCysCys'] = filterNumberCysCys
patterAnalysis['pattern'] = filterPattern

patterAnalysis.to_csv(sys.argv[2], index=False)
