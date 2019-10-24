'''
script que permite unir los elementos de la clasificacion con el fin de evaluar a la clase que corresponde
se generan diversos archivos que permiten
'''
import sys
import pandas as pd

#recibimos los dataFrame
clfH1 = pd.read_csv(sys.argv[1])
clfFull = pd.read_csv(sys.argv[2])
