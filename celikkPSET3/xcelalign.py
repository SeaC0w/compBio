# by kerim celik for computational biology
# 01/25/2018
import sys
import math
import numpy
import xlsxwriter
from Bio import SeqIO

# gets the index of the largest item in the passed list
# if multiple items are tied for the largest, return
# the first one in the list
def getMaxIndex(ls):
    return numpy.argmax(ls)

# cleanly prints a matrix
def printMatrix(m):
    for rows in m:
        print(rows)

# creates a matrix with the specified number of rows and columns, with
# all entries initialized to passed value
# uses local alignment algorithm to determine optimal alignment
def createMatrix(numRows, numColumns, value):
    matrix = []
    for a in range(numRows):
        ls = []
        for b in range(numColumns):
            ls.append(value)
        matrix.append(ls)
    return matrix

# takes in 2 strings v and w, and score list
# score list arranged as follows: match, mismatch, gap-open, gap-extend
# gap open score not used
def localAlign(v, w, scores):
    # empty string edge case handling
    if v == "" or w == "":
        print("Maximum Alignment Score: 0")
        print("--Alignment--")
        print("V")
        print("W")
        print("\n")
        return

    # attach gap to front of strings
    v = "-" + v
    w = "-" + w

    # create/initialize the dynamic programming table
    # also create a table to keep track of the traceback
    tabl = createMatrix(len(v), len(w), 0)
    trace = createMatrix(len(v), len(w), 3)
    # recursion step; keep track of max value in table
    maxEntry = (0,0)
    for i in range(1,len(v)):
        for j in range(1, len(w)):
            temp = []
            temp.append(tabl[i-1][j] + scores[3])
            temp.append(tabl[i][j-1] + scores[3])
            if v[i] == w[j]:
                temp.append(tabl[i-1][j-1] + scores[0])
            else:
                temp.append(tabl[i-1][j-1] + scores[1])
            temp.append(0)
            # get the max of the four localAlign recur values
            tabl[i][j] = max(temp)
            # traceback bookkeeping
            ind = getMaxIndex(temp)
            trace[i][j] = ind
            # keep track of the table's maximum value
            if temp[ind] > tabl[maxEntry[0]][maxEntry[1]]:
                maxEntry = (i,j)
    # traceback from maximum value in table to any table entry with a 0 in it
    maxAlignScore = tabl[maxEntry[0]][maxEntry[1]]
    currV = maxEntry[0] - 1
    currW = maxEntry[1] - 1
    retStrV = v[currV]
    retStrW = w[currW]
    # printMatrix(tabl)
    # print("\n")
    # printMatrix(trace)
    # while the current entry contains a value above 0
    workbook = xlsxwriter.Workbook("ratCatLocal.xlsx")
    worksheet = workbook.add_worksheet()
    while trace[currV][currW] != 3:
        worksheet.write(currV+1, currW+1, 1)
        # horizontal move
        if trace[currV][currW] == 1:
            retStrW = w[currW] + retStrW
            retStrV = "-" + retStrV
            currW -= 1
        # vertical move
        elif trace[currV][currW] == 0:
            retStrV = v[currV] + retStrV
            retStrW = "-" + retStrW
            currV -= 1
        # diagonal move
        else:
            retStrV = v[currV] + retStrV
            retStrW = w[currW] + retStrW
            currV -= 1
            currW -= 1
    worksheet.write(currV+1, currW+1, 1)
    worksheet.write(len(v)+1, len(w)+1,maxAlignScore)
    workbook.close()
    # return optimal alignment score and the two aligned strings
    print("Maximum Alignment Score: " + str(maxAlignScore))
    print("--Alignment--")
    print("V")
    print(retStrV)
    print("W")
    print(retStrW)

# takes in 2 strings v and w, and score list
# score list arranged as follows: match, mismatch, gap-open, gap-extend
# uses global alignment with affine gap algorithm to determine optimal alignment
def globalAlignAffineGap(v, w, scores):
    # empty string edge case handling
    if v == "" and w == "":
        print("Maximum Alignment Score: 0")
        print("--Alignment--")
        print("V")
        print("W")
        print("\n")
        return
    elif v == "":
        print("Maximum Alignment Score: " + str((scores[3] * len(w)) + scores[2]))
        retStr = ""
        for j in range(len(w)):
            retStr = retStr + "-"
        print("--Alignment--")
        print("V")
        print(retStr)
        print("W")
        print(w)
        return
    elif w == "":
        print("Maximum Alignment Score: " + str((scores[3] * len(v)) + scores[2]))
        retStr = ""
        for i in range(len(v)):
            retStr = retStr + "-"
        print("--Alignment--")
        print("V")
        print(v)
        print("W")
        print(retStr)
        return

    # attach gap to front of strings
    v = "-" + v
    w = "-" + w

    # create the dynamic programming tables
    # also create table to keep track of the traceback
    tablS = createMatrix(len(v), len(w), 0)
    tablE = createMatrix(len(v), len(w), 0)
    tablF = createMatrix(len(v), len(w), 0)
    tablG = createMatrix(len(v), len(w), 0)
    trace = createMatrix(len(v), len(w), -1)

    # initialize the tables, as specified by the algorithm
    for j in range(len(w)):
        tablE[0][j] = scores[2] + j * scores[3]
        tablS[0][j] = scores[2] + j * scores[3]
        tablF[0][j] = math.inf * -1
        tablG[0][j] = math.inf * -1
        trace[0][j] = 0
    for i in range(len(v)):
        tablF[i][0] = scores[2] + i * scores[3]
        tablS[i][0] = scores[2] + i * scores[3]
        tablE[i][0] = math.inf * -1
        tablG[i][0] = math.inf * -1
        trace[i][0] = 1
    tablS[0][0] = 0
    tablE[0][0] = scores[2]
    tablG[0][0] = 0
    trace[0][0] = 2

    # recursion step
    for i in range(1,len(v)):
        for j in range(1,len(w)):
            # calculate for table E
            temp = [tablE[i][j-1] + scores[3], tablS[i][j-1] + scores[2] + scores[3]]
            tablE[i][j] = max(temp)
            # calculate for table F
            temp = [tablF[i-1][j] + scores[3], tablS[i-1][j] + scores[2] + scores[3]]
            tablF[i][j] = max(temp)
            # calculate for table G
            if v[i] == w[j]:
                tablG[i][j] = tablS[i-1][j-1] + scores[0]
            else:
                tablG[i][j] = tablS[i-1][j-1] + scores[1]
            # calculate for table S
            temp = [tablE[i][j], tablF[i][j], tablG[i][j]]
            trace[i][j] = getMaxIndex(temp)
            tablS[i][j] = max(temp)

    # traceback from bottom-right corner to top-left corner
    currV = len(v) - 1
    currW = len(w) - 1
    retStrV = v[currV]
    retStrW = w[currW]
    workbook = xlsxwriter.Workbook("ratCatGlobal.xlsx")
    worksheet = workbook.add_worksheet()
    while currV != 0 or currW != 0:
        worksheet.write(currV+1, currW+1, 1)
        # horizontal move
        if trace[currV][currW] == 0:
            retStrW = w[currW] + retStrW
            retStrV = "-" + retStrV
            currW -= 1
        # vertical move
        elif trace[currV][currW] == 1:
            retStrV = v[currV] + retStrV
            retStrW = "-" + retStrW
            currV -= 1
        # diagonal move
        else:
            retStrV = v[currV] + retStrV
            retStrW = w[currW] + retStrW
            currV -= 1
            currW -= 1
    worksheet.write(currV+1, currW+1, 1)
    worksheet.write(len(v)+1, len(w)+1,tablS[len(v) - 1][len(w) - 1])
    workbook.close()
    # return optimal alignment score and the two aligned strings
    print("Maximum Alignment Score: " + str(tablS[len(v) - 1][len(w) - 1]))
    print("--Alignment--")
    print("V")
    print(retStrV)
    print("W")
    print(retStrW)

def main():
    # handles command line arguments, converting to usable forms
    fastaV = sys.argv[1]
    fastaW = sys.argv[2]
    scoreFile = sys.argv[3]
    option = int(sys.argv[4])
    seqV = SeqIO.parse(fastaV, "fasta")
    seqW = SeqIO.parse(fastaW, "fasta")
    v = ""
    w = ""
    for item in seqV:
        v = item.seq
    for item in seqW:
        w = item.seq
    v = v.lower()
    w = w.lower()
    print(len(v))
    f = open(scoreFile, 'r')
    fil = f.read()
    nums = fil.split("\n")[1]
    nums = nums.split("\t")
    delta = [int(item) for item in nums]
    # fixes my switching gap-open and gap-extend during implementation
    tempNum = delta[3]
    delta[3] = delta[2]
    delta[2] = tempNum
    # print(fil)
    print(delta)
    if option == 0:
        localAlign(v, w, delta)
    elif option == 1:
        globalAlignAffineGap(v, w, delta)
    f.close()

if __name__ == "__main__":
    main()
