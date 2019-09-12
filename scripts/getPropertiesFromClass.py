'''
script que permite obtener diferentes propiedades de las secuencias existentes en cada clase
genera un archivo JSON resumen con la informacion obtenida.

Permite obtener:

1. Largo maximo, minimo, promedio
2. Numero de cys maximo, minimo, promedio
3. Distribucion del largo del patron
4. Distribucion de residuos en el patron
'''
import sys
from Bio import SeqIO
import re

#obtenemos los parametros
fileData = sys.argv[1]
pathData = sys.argv[2]
typeClass = sys.argv[3]
