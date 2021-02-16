from constants import Matrix, shrunkMatrix
import cplex
import sys
import traceback
from itertools import compress
import time
import os



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


def isValidTrail(m):
    isValid = True
    SolLinearLayer = []
    for X, Y in m._InputOutputVariablesLinearLayer:
        sol = (m.solution.get_values(X), m.solution.get_values(Y))
        if sol not in SolLinearLayer:
            SolLinearLayer.append(sol)
    for u, v in SolLinearLayer:
        if sum(u + v) > 0:
            if len(u) != len(v):
                if not isNonSingular(subMatrix(shrunkMatrix, u, v)):
                    isValid = False
            else:
                if not isNonSingular(subMatrix(Matrix, u, v)):
                    isValid = False
    return isValid

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
            u, v = (list(map(round, context.get_candidate_point(X))), list(map(round, context.get_candidate_point(Y))))
            if len(X) != len(Y):
                if (sum(u + v) > 0) and (not isNonSingular(subMatrix(shrunkMatrix, u, v))):
                    c, s, r = self.excludeSolLazyCplex(X + Y, u + v)
                    _constraints.append(c)
                    _senses.append(s)
                    _rhs.append(r)
            else:
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



def SolveModelBitbyBitSingleBit(m):
    

    time_start = time.time()
    fileobj = open(m._filenameResult, "a", buffering=1)
    
    
    
    myCallBack = Callback(m._InputOutputVariablesLinearLayer)
    contextmask = cplex.callbacks.Context.id.candidate
    m.set_callback(myCallBack, contextmask)


    machine = os.uname()[1]
    set_zero = []
    for counter in range(0, 1):
        u = m._OutputVariables[counter]
        m.solve()
        status = m.solution.get_status()
        if status == m.solution.status.MIP_optimal:
            if isValidTrail(m):
                fileobj.write("%03i" %counter + " (" + u + ") " + " -- UNKNOWN\n")
                set_zero.append(u)
                m.write("sol/%03i_%s_%s.sol" %(m._shrunkColIndex, u, machine.split('.')[0]))
            else:
                fileobj.write("Lazy Constraints not work: Sol should be excluded. \n")
                m.write("Sol2Excluded_%i.sol" %counter)
                return([])

        elif status == m.solution.status.MIP_infeasible:
            fileobj.write("%03i" %counter + " (" + u + ") " + " -- Balanced\n")     
        else:
            fileobj.write("Unknown error!" + "\n")
            fileobj.close()
            return ([])        
    
    time_end = time.time()
    fileobj.close()

    balanceBits = []
    for i in range(len(m._OutputVariables)):
        varName = m._OutputVariables[i]
        if varName not in set_zero:
            balanceBits.append(varName)
    return (balanceBits)