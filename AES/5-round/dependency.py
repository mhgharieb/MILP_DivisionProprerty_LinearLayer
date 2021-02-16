from constants import shrunkMatrix
from itertools import combinations
from functools import reduce

def linearDependent(A):
    m = len(A)
    n = len(A[0])
    operationMat = [{i} for i in range(m)] 
    startRowIndx = 0
    for col in range(n):
        onesRowsIndx = [j for j in range(startRowIndx, m) if A[j][col] == 1]
        if len(onesRowsIndx) >= 1:
            pivotRowIndx = onesRowsIndx[0]
            pivotRow = A[pivotRowIndx]
            for c in range(1,len(onesRowsIndx)):
                currentRowIndx = onesRowsIndx[c]
                A[currentRowIndx] = [a ^ b for a, b in zip(A[currentRowIndx], pivotRow)]
                # print("*************************")
                # print(operationMat)
                # print(currentRowIndx, operationMat[currentRowIndx])
                # print(pivotRowIndx, operationMat[pivotRowIndx])
                # operationMat[currentRowIndx] = list(set(operationMat[currentRowIndx])^set(operationMat[pivotRowIndx]))
                operationMat[currentRowIndx] ^= operationMat[pivotRowIndx]
                # print(currentRowIndx, operationMat[currentRowIndx])
            
            # row swap
            tmp = list(pivotRow) 
            A[pivotRowIndx] = list(A[startRowIndx])
            A[startRowIndx] = list(tmp) 
            
            tmp = operationMat[pivotRowIndx]
            operationMat[pivotRowIndx] = operationMat[startRowIndx]
            operationMat[startRowIndx] = tmp 

            startRowIndx += 1
        
    DeprendentRowsSet_basis = []
    for rowIndx in range(m):
        if sum(A[rowIndx]) == 0:
            DeprendentRowsSet_basis.append(operationMat[rowIndx])
    
    DeprendentRowsSet_all = []
    NBasis = len(DeprendentRowsSet_basis)
    for i in range(1, NBasis + 1):
        for c in combinations(DeprendentRowsSet_basis, i):
            newDeprendentRowsSet = reduce(lambda x, y: x ^ y, c)
            toAppend = True
            for DeprendentRowsSet in DeprendentRowsSet_all:
                if DeprendentRowsSet.issubset(newDeprendentRowsSet):
                    toAppend = False
                    break
            if toAppend:
                DeprendentRowsSet_all.append(newDeprendentRowsSet)

    return DeprendentRowsSet_all

if __name__ == "__main__":
    a = linearDependent(shrunkMatrix)
    # print(len(a))
    print(a)
