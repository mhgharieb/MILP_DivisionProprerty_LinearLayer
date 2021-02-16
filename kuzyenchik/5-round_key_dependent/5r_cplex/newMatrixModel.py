# from constants import shrunkMatrixConstraints
from AESLike import words2bits


def ConstraintsByMixColumnNew(X_bytes, Y_bytes, r, ColIndex, MatrixConstraints):
    assert(len(X_bytes) == len(Y_bytes))
    X = words2bits(X_bytes)
    Y = words2bits(Y_bytes)
    Variables = []
    Constraints = []
    n = len(X)
    # sum(u) == sum(v)
    Constraints += [" + ".join(X) + ' - ' + " - ".join(Y) + " = 0"]
    
    Constraints +=[' + '.join([Y[i] for i in c]) + ' <= %i' %(len(c) - 1) for c in MatrixConstraints] 
    return (Variables, Constraints, (X, Y))    




