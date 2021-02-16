from constants import Matrix
import cplex
import sys
import traceback
from itertools import compress
import time
import os

from AESLike import words2bits, CreateStateVariable, MatrixTranspose, ShiftRow

def findAllUNKNOWN(m):
    Xout = CreateStateVariable('x', m._blocksize//m._wordSize, m._wordSize, m._Rounds)
    Z = ShiftRow(CreateStateVariable('y', m._blocksize//m._wordSize, m._wordSize, m._Rounds-1))
    Z_serialized = words2bits(Z)
    Z_UNKNOWN = Z_serialized[getVarsValue(m, Z_serialized).index(1)]


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
def subMatrix(M,u, v):
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



def getVarsValue(m, VarList):
    return list(map(round, m.solution.get_values(VarList)))


###################################################################
def isValidTrail(m):
    isValid = True
    SolLinearLayer = []
    for X, Y in m._InputOutputVariablesLinearLayer:
        sol = (getVarsValue(m, X), getVarsValue(m, Y))
        if sol not in SolLinearLayer:
            SolLinearLayer.append(sol)
    for u, v in SolLinearLayer:
        if sum(u + v) > 0:
            if not isNonSingular(subMatrix(Matrix, u, v)):
                isValid = False
    return isValid



###############################################################################
############Cplex CallBack#######################################

class Callback():
    def __init__(self, InputOutputVariablesLinearLayer):
        self.InputOutputVariablesLinearLayer = InputOutputVariablesLinearLayer

    def excludeSolLazyCplex(self, variables, mask):
        coef = [-2 * i + 1 for i in mask]
        senses = "G"
        rhs = 1 - sum(mask)
        return (cplex.SparsePair(variables, coef), senses, rhs)

    
    def isValidTrailCplex(self, context):
        
        if not context.is_candidate_point():
            raise Exception('Unbounded solution')
        
        _constraints = []
        _senses = []
        _rhs = []
        for X, Y in self.InputOutputVariablesLinearLayer:
            # m.cbGetSolution return float
            u, v = (list(map(round, context.get_candidate_point(X))), list(map(round, context.get_candidate_point(Y))))
            if (sum(u + v) > 0) and (not isNonSingular(subMatrix(Matrix, u, v))):
                c, s, r = self.excludeSolLazyCplex(X + Y, u + v)
                _constraints.append(c)
                _senses.append(s)
                _rhs.append(r)
        
        if len(_constraints) > 0:
            context.reject_candidate(
                    constraints=_constraints,
                    senses=_senses,
                    rhs=_rhs)

    def invoke(self, context):
        try:
            if context.in_candidate():
                self.isValidTrailCplex(context)
        except:
            info = sys.exc_info()
            print('#### Exception in callback: ', info[0])
            print('####                        ', info[1])
            print('####                        ', info[2])
            traceback.print_tb(info[2], file=sys.stdout)
            raise
#endif





def writeSol(m, SolFile):
    SolFile_obj = open(SolFile, "w")
    for name, value in zip(m.variables.get_names(),
                           m.solution.get_values()):
        SolFile_obj.write("%s %i\n" %(name, round(value)))
    SolFile_obj.close()



###########################################################
def SolveModelOptimize_findAllUNKNOWN(m):
    
    time_start = time.time()
    fileobj = open(m._filenameResult, "a", buffering=1)
    
    
    
    myCallBack = Callback(m._InputOutputVariablesLinearLayer)
    contextmask = cplex.callbacks.Context.id.candidate
    m.set_callback(myCallBack, contextmask)


    machine = os.uname()[1]
    set_zero = []
    set_zero_constraints = []
    counter = 0
    while (counter < m._blocksize):
        try:
            print("try: (CPX_MIPEMPHASIS_BALANCED, CPX_PARALLEL_OPPORTUNISTIC)", flush=True)
            m.parameters.emphasis.mip.set(0) #CPX_MIPEMPHASIS_BALANCED (default)
            m.parameters.parallel.set(-1) #CPX_PARALLEL_OPPORTUNISTIC   
            m.solve()
        except Exception as e:
            print("try: (CPX_MIPEMPHASIS_BALANCED, CPX_PARALLEL_OPPORTUNISTIC), error:%s" %(e), flush=True)
            try:
                print("try: (CPX_MIPEMPHASIS_BALANCED, CPX_PARALLEL_AUTO)", flush=True)
                m.parameters.emphasis.mip.set(0) #CPX_MIPEMPHASIS_FEASIBILITY
                m.parameters.parallel.set(0)
                m.solve()
            except Exception as e:
                print("try: (CPX_MIPEMPHASIS_BALANCED, CPX_PARALLEL_AUTO), error:%s" %(e), flush=True)
        
        status = m.solution.get_status()
        if status == m.solution.status.MIP_optimal:
            if isValidTrail(m):
                if round(m.solution.get_objective_value()) > 1:
                    break
                else:
                    u = m._OutputVariables[getVarsValue(m, m._OutputVariables).index(1)]
                    fileobj.write("%03i" %counter + " (" + u + ") " + " -- UNKNOWN\n")
                    print("%03i" %counter + " (" + u + ") " + " -- UNKNOWN\n")
                    writeSol(m, "sol/%s_%s_%s.sol" %(m._SolNamePrefix, u, machine.split('.')[0]))

                    Xout_AllUNKNOWN = findAllUNKNOWN(m)
                    toSetZero = [v for v in  Xout_AllUNKNOWN if v not in set_zero]
                    set_zero += list(toSetZero)

                    for u in toSetZero:
                        fileobj.write("%03i" %counter + " (" + u + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(0))
                        print("%03i" %counter + " (" + u + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(0))
                        set_zero_constraints += m.linear_constraints.add(
                            lin_expr = [cplex.SparsePair(ind = [u], val = [1.0])],
                            senses = ["E"],
                            rhs = [0.0],
                            names = ['tmp%i'%counter])
                        
                        counter += 1
            else:
                fileobj.write("Lazy Constraints not work: Sol should be excluded. \n")
                writeSol(m, "Sol2Excluded_%s_%s_%s.sol" %(m._SolNamePrefix, u, machine.split('.')[0]))
                return([])

        elif status == m.solution.status.MIP_infeasible:
            break
        else:
            fileobj.write("Unknown error!" + "\n")
            # print("Unknown error!" + "\n")
            fileobj.close()
            return ([])        
    
    time_end = time.time()
    fileobj.write("run time: %f Seconds \n" %(time_end - time_start) )
    fileobj.close()


    if len(set_zero_constraints) > 0:
        m.linear_constraints.delete(set_zero_constraints)


    balanceBits = []
    for i in range(len(m._OutputVariables)):
        varName = m._OutputVariables[i]
        if varName not in set_zero:
            balanceBits.append(varName)
    return (balanceBits)    