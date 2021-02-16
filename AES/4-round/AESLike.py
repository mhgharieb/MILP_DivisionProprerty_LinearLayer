from constants import *


def CreateStateVariable(x, N, wordSize, Track):
    state = []
    for iByte in range(N):
        state.append([x + '_%02i' % Track + '_%02i' %
                      iByte + '_%i' % j for j in reversed(range(wordSize))])
    return state       

def CreateTmpVariable(x, s):
    return [s+ "_%03i_%03i" %(r,c) for r,c in x]


def words2bits(X):
    return sum(X, [])  

def Copy(a, B):
    #(a) --COpy--> (b0, b1, ... , bm)
    return [a + ' - ' + ' - '.join(B) + ' = 0']

def XOR(A, b):
    #(a0, a1, ... , am) --XOR--> (b)
    return [' + '.join(A) + ' - ' + b + ' = 0']


def CopyArray(A):
    copyArray = []
    for i in range(len(A)):
        tmp = [(j,i) for j, x in enumerate(A[i]) if x == 1]
        copyArray.append(tmp)
    return copyArray

def XorArray(A):
    xorArray = []
    for i in range(len(A)):
        tmp = [(i,j) for j, x in enumerate(A[i]) if x == 1]
        xorArray.append(tmp)    
    return xorArray


# Y = LinearMatrix * X , where X = [x0, x1, ..., xn] and x0 is LSB 
MatrixTranspose = list(map(list, zip(*Matrix)))

CopyMatrix =  CopyArray(MatrixTranspose)
XorMatrix = XorArray(Matrix)


shrunkMatrixTranspose = list(map(list, zip(*shrunkMatrix)))

CopyshrunkMatrix =  CopyArray(shrunkMatrixTranspose)
XorshrunkMatrix = XorArray(shrunkMatrix)


def ConstraintsByMixColumn(X_bytes, Y_bytes, r, ColIndex):
    # print(len(X_bytes))
    # print(len(Y_bytes))
    assert(len(X_bytes) == len(Y_bytes))
    X = words2bits(X_bytes)
    Y = words2bits(Y_bytes)
    Variables = []
    Constraints = []
    n = len(X)
    # sum(u) == sum(v)
    Constraints += [" + ".join(X) + ' - ' + " - ".join(Y) + " = 0"]
    #copyConstraints
    for vPos in range(n):
        V = CreateTmpVariable(CopyMatrix[vPos], 't_%02i_%02i' %(r, ColIndex))
        Variables += V
        Constraints += Copy(X[vPos], V)
    #xorConstraints
    for vPos in range(n):
        V = CreateTmpVariable(XorMatrix[vPos], 't_%02i_%02i' %(r, ColIndex))
        #Variables += V
        Constraints += XOR(V, Y[vPos])
    return (Variables, Constraints, (X, Y))

def ConstraintsByShrunkMixColumn(X_bytes, Y_bytes, r, ColIndex):
    assert(len(X_bytes) + 1 == len(Y_bytes))
    X = words2bits(X_bytes)
    Y = words2bits(Y_bytes)
    Variables = []
    Constraints = []
    n = len(X)
    # sum(u) == sum(v)
    Constraints += [" + ".join(X) + ' - ' + " - ".join(Y) + " = 0"]
    #copyConstraints
    for vPos in range(n):
        V = CreateTmpVariable(CopyshrunkMatrix[vPos], 't_%02i_%02i' %(r, ColIndex))
        Variables += V
        Constraints += Copy(X[vPos], V)
    #xorConstraints
    for vPos in range(n):
        V = CreateTmpVariable(XorshrunkMatrix[vPos], 't_%02i_%02i' %(r, ColIndex))
        #Variables += V
        Constraints += XOR(V, Y[vPos])
    return (Variables, Constraints, (X, Y))

def ConstraintsByLinearLayer(X, Y, ColSize, r):
    Variables = []
    Constraints = []
    InputOutputVariable = []
    for ColIndex in range(len(X)//ColSize):
        XCurrentCol = X[ColSize*ColIndex:ColSize*(ColIndex+1)]
        YCurrentCol = Y[ColSize*ColIndex:ColSize*(ColIndex+1)]
        variables , constraints, tmpInputOutputVariable = ConstraintsByMixColumn(XCurrentCol, YCurrentCol, r, ColIndex)
        Variables += variables
        Constraints +=  constraints
        InputOutputVariable.append(tmpInputOutputVariable)
    return (Variables, Constraints, InputOutputVariable)


def ConstraintsBySbox(X, Y):
    V = X + Y
    constraints = []
    for inq in SboxInqs:
        tmp = [str(a) + ' ' + str(b) for (a, b) in zip(inq[:-1], V)]
        tmp1 = ' + '.join(tmp)
        tmp2 = tmp1.replace("+ -", "- ")
        tmp2 += ' >= '
        tmp2 += str(-1 * inq[-1])
        constraints.append(tmp2)

    return constraints
    

def ConstraintsBySboxLayer(Xbytes, Ybytes):
    constraints = []
    n = len(Xbytes)
    for i in range(n):
        constraints += ConstraintsBySbox(Xbytes[i], Ybytes[i])
    return constraints


def ShiftRow(X):
    return [X[i] for i in ShiftRowIndex]


 
    
