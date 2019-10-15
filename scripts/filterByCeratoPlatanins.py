'''
script que permite aplicar el filtro de motivos conservados en las secuencias para proteinas del tipo
cerato platanins

Los filtros aplicados son:

1. Numero de cys
2. Numero de cys-cys
3. Patron caracteristico clase II
'''

import sys
from Bio import SeqIO
import pandas as pd
import re

#funcion que permite contar el numero de cys, debe ser mayor o igual a 4 y menor a 6
def countNumberCys(sequence):

    contCys=0

    for i in range(len(sequence)):
        if sequence[i] == 'C':
            contCys+=1
    if contCys>=4 and contCys<=6:#cumple con el minimo
        return 0
    else:#no lo cumple
        return 1

#funcion que permite contar los cys-cys, debe ser mayor o igual a dos
def countNumberCysCys(sequence):

    contCysCys = 0

    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            contCysCys+=1
            break

    return contCysCys#ya que no debe tener dobles cys

#funcion que permite evaluar el patron de clase I reportado en la literatura
def evaluatedPattern(sequence, patternString):

    pattern = re.compile(patternString)
    response = pattern.findall(str(sequence))

    if len(response)>0:#cumple con el patron
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

#seccion de patrones
filterPattern1 = []
contFilterPattern1 = 0

filterPattern2 = []
contFilterPattern2 = 0

filterPattern3 = []
contFilterPattern3 = 0

filterPattern4 = []
contFilterPattern4 = 0

filterPattern5 = []
contFilterPattern5 = 0

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

    #analizamos patron VS
    response3 = evaluatedPattern(record.seq, 'VS')
    filterPattern1.append(response3)
    if response3 == 0:
        contFilterPattern1+=1

    #analizamos patron YD[XXXXX]YD
    response4 = evaluatedPattern(record.seq, 'YD[A-Z]{3,7}YD')
    filterPattern2.append(response4)
    if response4 == 0:
        contFilterPattern2+=1

    #analizamos patron C[XX]S[X]D[XX]G
    response5 = evaluatedPattern(record.seq, 'C[A-Z]{1,3}S[A-Z]{1,2}D[A-Z]{1,3}G')
    filterPattern3.append(response5)
    if response5 == 0:
        contFilterPattern3+=1

    #analizamos patron C[XX]G
    response6 = evaluatedPattern(record.seq, 'C[A-Z]{1,3}G')
    filterPattern4.append(response6)
    if response6 == 0:
        contFilterPattern4+=1

    #analizamos patron VLAID
    response7 = evaluatedPattern(record.seq, 'VLAID')
    filterPattern5.append(response7)
    if response7 == 0:
        contFilterPattern5+=1

print cont_filterNumberCys
print cont_filterNumberCysCys
print contFilterPattern1
print contFilterPattern2
print contFilterPattern3
print contFilterPattern4
print contFilterPattern5

#agregamos la data al dataFrame y exportamos los analisis
patterAnalysis['numberCys'] = filterNumberCys
patterAnalysis['numberCysCys'] = filterNumberCysCys
patterAnalysis['pattern1'] = filterPattern1
patterAnalysis['pattern2'] = filterPattern2
patterAnalysis['pattern3'] = filterPattern3
patterAnalysis['pattern4'] = filterPattern4
patterAnalysis['pattern5'] = filterPattern5

patterAnalysis.to_csv(sys.argv[2], index=False)
