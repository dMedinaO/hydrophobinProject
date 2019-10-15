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
    if contCys>=8:#cumple con el minimo
        return 0
    else:#no lo cumple
        return 1

#funcion que permite contar los cys-cys, debe ser mayor o igual a dos
def countNumberCysCys(sequence):

    contCysCys = 1

    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            contCysCys=0
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
    response3 = evaluatedPattern(record.seq, 'PQCC')
    filterPattern1.append(response3)
    if response3 == 0:
        contFilterPattern1+=1

print cont_filterNumberCys
print cont_filterNumberCysCys
print contFilterPattern1

#agregamos la data al dataFrame y exportamos los analisis
patterAnalysis['numberCys'] = filterNumberCys
patterAnalysis['numberCysCys'] = filterNumberCysCys
patterAnalysis['pattern1'] = filterPattern1

patterAnalysis.to_csv(sys.argv[2], index=False)
