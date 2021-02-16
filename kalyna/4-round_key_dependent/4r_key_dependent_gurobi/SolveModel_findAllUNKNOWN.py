from constants import Matrix
from gurobipy import GRB
from itertools import compress
import time
import os



from AESLike import words2bits, CreateStateVariable, MatrixTranspose, ShiftRow

def findAllUNKNOWN(m):
    Xout = CreateStateVariable('x', m._blocksize//m._wordSize, m._wordSize, m._Rounds)
    Z = ShiftRow(CreateStateVariable('y', m._blocksize//m._wordSize, m._wordSize, m._Rounds-1))
    Z_serialized = words2bits(Z)
    Z_UNKNOWN = Z_serialized[getVarsValue(getVarsByName(m, Z_serialized)).index(1)]


    Z_column = [Z[i:i+m._ColSize] for i in range(0, m._blocksize//m._wordSize, m._ColSize)]
    Xout_column = [Xout[i:i+m._ColSize] for i in range(0, m._blocksize//m._wordSize, m._ColSize)]
    for ColIndex in range(len(Z_column)):
        currentZcol = words2bits(Z_column[ColIndex])
        currentXoutcol = words2bits(Xout_column[ColIndex]) 
        if Z_UNKNOWN in currentZcol:
            Z_UNKNOWN_Indx = currentZcol.index(Z_UNKNOWN)
            allUNKNOWN = list(compress(currentXoutcol, MatrixTranspose[Z_UNKNOWN_Indx])) # column in Matrix
            return allUNKNOWN



#--------------------------------------------------------------------------------------
def subMatrix(M, u, v):
    return [list(compress(row, u)) for row, v_item in zip(M, v) if v_item == 1]

def isNonSingular(A):
    n = len(A[0])
    m = n
    for col in range(n):
        rows = [j for j in range(m) if A[j][col] == 1]
        if len(rows) >= 1:
            row = A[rows[0]]
            for c in range(1,len(rows)):
                A[rows[c]] = [a ^ b for a, b in zip(A[rows[c]], row)]
            del A[rows[0]]
            m -= 1
        else:
            return False
    return True

def getVarsByName(m, VarList):
    return [m.getVarByName(var) for var in VarList]    

def getVarsValue(VarList):
    return [round(var.x) for var in VarList]

def setInputDP(m, inputDP):
    inputConstraints = m.addConstrs((m._InputVariables[i] == int(inputDP[i]) for i in range(len(inputDP))), name='inputConstraints')
    m.update()
    return inputConstraints

def resetOutputBound(m):
    for u in m._OutputVariables:
        u.ub = 1
        u.lb = 0
    m.update()
###################################################################
def isValidTrail(m):
    isValid = True
    SolLinearLayer = []
    for X, Y in m._InputOutputVariablesLinearLayer:
        sol = (getVarsValue(X), getVarsValue(Y))
        if sol not in SolLinearLayer:
            SolLinearLayer.append(sol)
    for u, v in SolLinearLayer:
        if sum(u + v) > 0:
            if sum(u) != sum(v):
                isValid = False
            elif not isNonSingular(subMatrix(Matrix, u, v)):
                isValid = False
    return isValid

def excludeSol(m, mask):
    print("invalid: " + str(mask), flush=True)
    OnesMsk = mask
    ZerosMsk = [not i for i in mask]
    for X, Y in m._InputOutputVariablesLinearLayer:
        variables = X + Y
        m.addConstr(sum(list(compress(variables, ZerosMsk))) + sum(OnesMsk) - sum(list(compress(variables, OnesMsk))) >= 1)
    m.update()



def getVarsName(X):
    return [u.getAttr('VarName') for u in X]

def excludeSolLazy(m, variables, mask):
    # print("invalid: " + str(mask) + ", " + str(getVarsName(variables)), flush=True)
    OnesMsk = mask
    ZerosMsk = [not i for i in mask]
    m.cbLazy(sum(list(compress(variables, ZerosMsk))) + sum(OnesMsk) - sum(list(compress(variables, OnesMsk))) >= 1)


def isValidTrailCallback(m, where):
    if where == GRB.Callback.MIPSOL:
        for X, Y in m._InputOutputVariablesLinearLayer:
            u, v = (list(map(round, m.cbGetSolution(X))), list(map(round, m.cbGetSolution(Y))))
            if (sum(u + v) > 0) and (not isNonSingular(subMatrix(Matrix, u, v))):
                excludeSolLazy(m, X + Y, u + v)


########################################################################################
def SolveModelOptimize_findAllUNKNOWN(m):
    

    time_start = time.time()
    fileobj = open(m._filenameResult, "a", buffering=1)
    
    m.setParam("LazyConstraints", 1)
    set_zero = []
    set_zero_constraints = []
    counter = 0
    while (counter < m._blocksize):
        m.optimize(isValidTrailCallback)
        # Gurobi syntax: m.Status == 2 represents the model is feasible.
        if m.Status == 2:
            if isValidTrail(m):
                obj = m.getObjective()
                if round(obj.getValue()) > 1:
                    break
                else:
                    for i in range(0, m._blocksize):
                        u = obj.getVar(i)
                        temp = round(u.getAttr('x'))
                        if temp == 1:
                            fileobj.write("%03i" %counter + " (" + u.getAttr('VarName') + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(m.Runtime))
                            print("%03i" %counter + " (" + u.getAttr('VarName') + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(m.Runtime))
                            m.write("sol/%s.sol" %u.getAttr('VarName'))
                            Xout_AllUNKNOWN = findAllUNKNOWN(m)
                            toSetZero = [v for v in  Xout_AllUNKNOWN if v not in set_zero]
                            set_zero += toSetZero

                            for u in getVarsByName(m, toSetZero):
                                fileobj.write("%03i" %counter + " (" + u.getAttr('VarName') + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(0))
                                print("%03i" %counter + " (" + u.getAttr('VarName') + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(0))
                                set_zero_constraints.append(m.addConstr( u == 0, 'tmp%i'%counter))
                                counter += 1
                            m.update()
                            break
            else:
                fileobj.write("Lazy Constraints not work: Sol should be excluded. \n")
                print("Lazy Constraints not work: Sol should be excluded. \n")
                m.write("Sol2Excluded_%i.sol" %counter)
                return([])

        # Gurobi syntax: m.Status == 3 represents the model is infeasible.
        elif m.Status == 3:
            break
        # Gurobi syntax: m.Status == 9 represents the model is TimeLimit.
        elif m.Status == 9:
            fileobj.write("%03i" %counter + " (" + u.getAttr('VarName') + ") " + " -- TimeLimit \n")
            print("%03i" %counter + " (" + u.getAttr('VarName') + ") " + " -- TimeLimit \n")
            set_zero.append(u.getAttr('VarName'))      
        else:
            fileobj.write("Unknown error!" + "\n")
            print("Unknown error!" + "\n")
            fileobj.close()
            return ([])        
    
    time_end = time.time()
    fileobj.write("run time: %f Seconds \n" %(time_end - time_start) )
    fileobj.close()
    
    if len(set_zero_constraints) > 0:
        m.remove(set_zero_constraints)
        m.update()

    balanceBits = []
    for i in range(len(m._OutputVariables)):
        varName = m._OutputVariables[i].getAttr('VarName')
        if varName not in set_zero:
            balanceBits.append(varName)
    return (balanceBits)        