from constants import Matrix
from gurobipy import GRB
from itertools import compress
import time
import os



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
    OnesMsk = mask
    ZerosMsk = [not i for i in mask]
    m.cbLazy(sum(list(compress(variables, ZerosMsk))) + sum(OnesMsk) - sum(list(compress(variables, OnesMsk))) >= 1)


def isValidTrailCallback(m, where):
    if where == GRB.Callback.MIPSOL:
        # SolLinearLayer = []
        for X, Y in m._InputOutputVariablesLinearLayer:
            # m.cbGetSolution return float
            u, v = (list(map(round, m.cbGetSolution(X))), list(map(round, m.cbGetSolution(Y))))
            if (sum(u + v) > 0) and (not isNonSingular(subMatrix(Matrix, u, v))):
                excludeSolLazy(m, X + Y, u + v)

#########################################################################
def SolveModelOptimize(m):

    time_start = time.time()
    fileobj = open(m._filenameResult, "a", buffering=1)
    
    m.setParam("LazyConstraints", 1)
    set_zero = []
    set_zero_constraints = []
    for counter in range(0, m._blocksize):
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
                            set_zero.append(u.getAttr('VarName'))
                            set_zero_constraints.append(m.addConstr(u == 0, 'tmp%i'%counter))
                            m.write("sol/%s.sol" %u.getAttr('VarName'))
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


