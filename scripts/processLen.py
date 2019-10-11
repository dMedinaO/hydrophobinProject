import pandas as pd
import sys
from Bio import SeqIO

inputFile = sys.argv[1]
pathOutput = sys.argv[2]

lenElements = []

for record in SeqIO.parse(inputFile, "fasta"):

    lenElements.append(len(record.seq))


dataFrame = pd.DataFrame(lenElements, columns=['lenData'])
dataFrame.to_csv(pathOutput+"lenData.csv", index=False)
