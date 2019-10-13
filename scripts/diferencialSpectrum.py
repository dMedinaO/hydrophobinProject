'''
script que recibe dos archivos de espectros y permite calcular los espectros diferenciales entre ellos
'''

import sys
import pandas as pd

#hacemos la diferencia de los espectros
def diffSpectrum(spectra1, spectra2):

    diffArray = []

def readDocument(document):

    docSpectrum = []
    fileDoc = open(document, 'r')
    line = fileDoc.readline()

    while line:
        row = line.replace("\n", "").split(",")
        docSpectrum.append(row)
        line = fileDoc.readline()
    fileDoc.close()

    return docSpectrum

docSpectrum1 = sys.argv[1]
docSpectrum2 = sys.argv[2]

spectrumData1 = readDocument(docSpectrum1)
spectrumData2 = readDocument(docSpectrum2)

print len(spectrumData1[0])
print len(spectrumData1[1])
