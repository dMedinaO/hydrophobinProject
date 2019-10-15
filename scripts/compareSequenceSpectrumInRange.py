'''
script que recibe el conjunto de espectros obtendos a traves de la digitalizacion de las propiedades fisicoquimicas de las secuencias
y ademas el archivo index con la data de interes, junto con el archivo de valores caracteristicos de los espectros de las secuencias
ya categorizadas y permite comparar los valores por medio de la evaluacion de zonas y obtiene porcentajes de los elementos que pertenecen a cada zona

Z1: Punto es mayor al espectro maximo
Z2: Punto mayor a IC max pero menor o igual a espectro maximo
Z3: Punto menor o igual a IC max y mayor o igual a IC min ------>>>>>>> de preferencia un porcentaje importante debiese estar en esta zona
Z4: Punto menor a IC min pero mayor o igual a espectro min
Z5: Punto menor a espectro min

Luego se evaluan las razones de los elementos y se reportan

Se recibe el archivo de procesamiento asociado a los IDs de las secuencias evaluadas, ademas de la propiedad a ser evaluada

Archivo de propiedades:
Pos[0] => min
Pos[1] => max
Pos[2] => mean
Pos[3] => IC min
Pos[4] => IC max
'''

import sys
import pandas as pd

#funcion que permite contar la cantidad de tipos de zona y retorna un arreglo con el contador
def countElementZone(zoneArray):

    countZ1 = 0
    countZ2 = 0
    countZ3 = 0
    countZ4 = 0
    countZ5 = 0

    for zone in zoneArray:
        if zone == 'Z1':
            countZ1+=1
        if zone == 'Z2':
            countZ2+=1
        if zone == 'Z3':
            countZ3+=1
        if zone == 'Z4':
            countZ4+=1
        if zone == 'Z5':
            countZ5+=1

    countZ1 = float(countZ1)*100/len(zoneArray)
    countZ2 = float(countZ2)*100/len(zoneArray)
    countZ3 = float(countZ3)*100/len(zoneArray)
    countZ4 = float(countZ4)*100/len(zoneArray)
    countZ5 = float(countZ5)*100/len(zoneArray)

    return countZ1, countZ2, countZ3, countZ4, countZ5

#funcion que permite asignar una zona a un punto
def defineZoneByPoint(pointSpectrum, minPoint, maxPoint, icMinPoint, icMaxPoint):

    zone = ''

    if pointSpectrum> maxPoint:
        zone='Z1'
    elif pointSpectrum<= maxPoint and pointSpectrum > icMaxPoint:
        zone='Z2'
    elif pointSpectrum <= icMaxPoint and pointSpectrum>= icMinPoint:
        zone='Z3'
    elif pointSpectrum <icMinPoint and pointSpectrum >= minPoint:
        zone='Z4'
    else:
        zone='Z5'

    return zone

#funcion que permite castear la data de string to float
def castElementsArray(arrayData):

    arrayCast = []

    for element in arrayData:
        arrayCast.append(float(element))
    return arrayCast

#funcion que permite la lectura y transformacion de la data... genera una matriz con los elementos
def readFileData (dataToRead):

    matrixResponse = []

    fileOpen = open(dataToRead, 'r')
    line = fileOpen.readline()

    while line:
        row = line.replace("\n", "").split(",")
        matrixResponse.append(castElementsArray(row))
        line = fileOpen.readline()
    return matrixResponse

#recibimos la data de interes
spectrumData = sys.argv[1]
idSequence = pd.read_csv(sys.argv[2])
spectrumProperties = sys.argv[3]
property = sys.argv[4]

#leemos y procesamos la data de los espectros
spectrumDataMatrix = readFileData(spectrumData)
spectrumPropertiesMatrix = readFileData(spectrumProperties)

#definimos arrays para agregar la proporcion de las zonas
zone1Array = []
zone2Array = []
zone3Array = []
zone4Array = []
zone5Array = []

#trabajamos con todos los elementos
for i in range (1, len(spectrumDataMatrix)):
    rowZone = []

    for j in range(len(spectrumDataMatrix[i])):#trabajamos con el espectro i con todos sus puntos
        rowZone.append(defineZoneByPoint(spectrumDataMatrix[i][j], spectrumPropertiesMatrix[0][j], spectrumPropertiesMatrix[1][j], spectrumPropertiesMatrix[3][j], spectrumPropertiesMatrix[4][j]))

    #contamos los elementos y la proprocion correspondiente
    countZ1, countZ2, countZ3, countZ4, countZ5 = countElementZone(rowZone)

    zone1Array.append(countZ1)
    zone2Array.append(countZ2)
    zone3Array.append(countZ3)
    zone4Array.append(countZ4)
    zone5Array.append(countZ5)

#formamos los nuevos elementos del data frame asociados con la key de interes
idSequence[property+"Z1"] = zone1Array
idSequence[property+"Z2"] = zone2Array
idSequence[property+"Z3"] = zone3Array
idSequence[property+"Z4"] = zone4Array
idSequence[property+"Z5"] = zone5Array

#exportamos el dataframe
idSequence.to_csv(sys.argv[2], index=False)
