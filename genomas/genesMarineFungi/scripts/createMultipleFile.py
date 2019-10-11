'''
script que recibe un archivo multifasta y crea archivos multifasta de menor cantidad de secuencias, con el fin de poder
ser procesados en alineamientos contra blast
'''

import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW

idsData = []
seqData = []

for record in SeqIO.parse(sys.argv[1], "fasta"):

    result_handle = NCBIWWW.qblast("blastp", "nr", record.seq, format_type="Text")
    with open(record.id, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
    break
