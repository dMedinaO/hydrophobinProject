import sys
from Bio import SeqIO

nameDoc = sys.argv[1]
posInit = int(sys.argv[2])
posFin = int(sys.argv[3])
idSegment = sys.argv[4]

sequenceData = ""
for record in SeqIO.parse(nameDoc, "fasta"):

    if record.id == idSegment:
        sequenceData = record.seq
        break
sequenceInGenome = ""

for i in range(posInit-1, posFin+1):

    sequenceInGenome+= sequenceData[i]

print sequenceInGenome
