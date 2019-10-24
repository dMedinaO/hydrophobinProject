'''
script que permite unir los elementos de la clasificacion con el fin de evaluar a la clase que corresponde
se generan diversos archivos que permiten
'''
import sys
import pandas as pd

#recibimos los dataFrame
clfH1 = pd.read_csv(sys.argv[1])
clfFull = pd.read_csv(sys.argv[2])
pathOutout = sys.argv[3]

#obtengo todos los ID de las secuencias
idSequences = []

for element in clfFull['indexSequence']:
    idSequences.append(str(element))
for element in clfH1['indexSequence']:
    idSequences.append(str(element))

idSequences = list(set(idSequences))

arrayClassCeratoP = []
arrayClassCeratoU = []
arrayClassH2 = []
arrayClassH1 = []

for idElement in idSequences:

    cont = 0
    #buscamos las clases en el FULL
    for i in range(len(clfFull)):
        if clfFull['indexSequence'][i] == idElement:
            arrayClassCeratoP.append(clfFull['ceratoPlatanin'][i])
            arrayClassCeratoU.append(clfFull['ceratoUlmin'][i])
            arrayClassH2.append(clfFull['hidroCII'][i])
            cont=1
            break
    #significa que no esta en ninguna clase
    if cont==0:
        arrayClassCeratoP.append(1)
        arrayClassCeratoU.append(1)
        arrayClassH2.append(1)

    cont2=0

    #buscamos la clase en la H1
    for i in range(len(clfH1)):
        if clfH1['indexSequence'][i] == idElement:
            arrayClassH1.append(clfH1['class'][i])
            cont2=1
            break
    #significa que no esta en ninguna clase
    if cont2==0:
        arrayClassH1.append(1)

#formamos el dataFrame
dataFull = pd.DataFrame()
dataFull['indexSequence'] = idSequences
dataFull['classCeratoP'] = arrayClassCeratoP
dataFull['classCeratoU'] = arrayClassCeratoU
dataFull['classH2'] = arrayClassH2
dataFull['classH1'] = arrayClassH1

#exportamos el dataFrame
dataFull.to_csv(pathOutout+"FullClassData.csv", index=False)
