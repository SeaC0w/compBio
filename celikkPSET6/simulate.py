# blurb
import sys
import math
import numpy
from Bio import SeqIO
import random

def simulateReads(readsFile, coverage, readLen, errorRate):
    readLen = int(readLen)
    coverage = int(coverage)
    errorRate = float(errorRate)
    rawSeq = SeqIO.parse(readsFile, 'fasta')
    seqToSim = ''
    for item in rawSeq:
        seqToSim = item.seq
    # print(seqToSim)
    seqToSim = seqToSim.lower().upper()

    seqToSimLen = len(seqToSim)
    # print("sequence length = " +  str(seqToSimLen))
    # print("coverage = " + str(coverage
    # print()
    if readLen > seqToSimLen:
        print("ERROR: your read length is longer than the sequence you want to read from.")
        return
    if ((0.0 >= errorRate) or (1.0 <= errorRate)):
        print("ERROR: your error rate was not between 0 and 1.")
    numReads = 0
    maxReads = (coverage * seqToSimLen) / readLen
    print(coverage)
    print(seqToSimLen)
    print(readLen)
    print("maxReads = " + str(maxReads))
    readsToWrite = []
    while numReads < maxReads:
        startIndex = random.randint(0, seqToSimLen - readLen + 1)
        currRead = seqToSim[startIndex:startIndex + readLen + 1]
        readsToWrite.append(currRead)
        numReads += 1
    # print(readsToWrite)
    fil = open("simulatedReads.txt", "w")
    readsFinal = [str(i) + '\n' for i in readsToWrite]
    fil.writelines(readsFinal)
    fil.close()


def main():
    if len(sys.argv) != 5:
        print("Error, invalid number of arguments submitted.")
        print("Should be: scriptName, fastaFilename, coverageInteger, readLengthInteger, errorRateFloat")
    rList = sys.argv[1]
    cov = sys.argv[2]
    rLen = sys.argv[3]
    eRate = sys.argv[4]
    simulateReads(rList, cov, rLen, eRate)

if __name__ == "__main__":
    main()
