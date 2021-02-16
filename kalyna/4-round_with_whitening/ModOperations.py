#! /usr/bin/env python


def CopyVariable(x, l):
    return [x + '_%02i' %i for i in range(l)]

def Copy(a, B):
    #(a) --COpy--> (b0, b1, ... , bm)
    return [a + ' - ' + ' - '.join(B) + ' = 0']

def XOR(A, b):
    #(a0, a1, ... , am) --XOR--> (b)
    return [' + '.join(A) + ' - ' + b + ' = 0']

def AND(A, b):
    #(a0, a1) --AND--> (b)
    return [
    b + ' - ' + A[0] + ' >= 0',
    b + ' - ' + A[1] + ' >= 0',
    b + ' - ' + A[0] + ' - ' + A[1] + ' <= 0'
    ]


def CreateVariables(x, Track, l):
    return [x + '_%02i' %Track + '_%02i' %i for i in range(l)]

def XORlist(A, B, C):
    #(A, B) --XOR--> (C)
    #([a0, a1, ... , am], [b0, b1, ... , bm]) --XOR--> ([c0, c1, ... , cm])
    #(a0, b0) --XOR--> (c0) , (a1, b1) --XOR--> (c1), ... , (am, bm) --XOR--> (cm)
    constraints = []
    for a, b, c in zip(A, B, C):
        constraints += XOR([a, b], c)
    return constraints

#def RotateLeft(X, m):
#    return X[m:] + X[:m]

#def RotateRight(X, m):
#    return X[-1*m:] + X[:-1*m]


def Copylist(A, B, C):
    #(A) --Copy--> (B, C)
    #([a0, a1, ... , am]) --Copy--> ()[b0, b1, ... , bm], [c0, c1, ... , cm])
    #(a0) --Copy--> (b0, c0) , (a1) --Copy--> (b1, c1), ... , (am) --Copy--> (bm, cm)
    constraints = []
    for a, b, c in zip(A, B, C):
        constraints += Copy(a, [b, c])
    return constraints




def ConstraintsByModAddConstant(X, Z, addTrack):
    assert(len(X) == len(Z))
    Variables = []
    Constraints = []
    n = len(X)
    #auxiliary variables
    v = CreateVariables('vA', addTrack, n-2)
    g = CreateVariables('gA', addTrack, n-2)
    f = CreateVariables('fA', addTrack, n-2)
    e = CreateVariables('eA', addTrack, n-2)

    Variables += v + g + f + e

    #1. (x_n-1) --Copy--> (z_n-1,f_0, g_0)
    Constraints += Copy(X[n-1], [Z[n-1], f[0], g[0]])

    #2. (x_n-2) --Copy--> (x_n-2_0, x_n-2_1, x_n-2_2)
    [a0, a1, a2] = CopyVariable(X[n-2], 3)
    Variables += [a0, a1, a2]
    Constraints += Copy(X[n-2], [a0, a1, a2])

    #3. (x_n-2_0, f_0) --XOR--> (z_n-2)
    Constraints += XOR([a0, f[0]], Z[n-2])

    #4. (x_n-2_1, g_0) --AND--> (e_0)
    Constraints += AND([a1, g[0]], e[0])

    #5. (x_n-2_2, e_0) --XOR--> (v_0)
    Constraints += XOR([a2, e[0]], v[0])

    for i in range(1, n-2):

        #6. (v_i-1) --Copy--> (f_i, g_i)
        Constraints += Copy(v[i-1], [f[i], g[i]])

        #7. (x_n-2-i) --Copy--> (x_n-2-i_0, x_n-2-i_1, x_n-2-i_2)
        [a0, a1, a2] = CopyVariable(X[n-2-i], 3)
        Variables += [a0, a1, a2]
        Constraints += Copy(X[n-2-i], [a0, a1, a2])

        #8. (x_n-2-i_0, f_i) --XOR--> (z_n-2-i)
        Constraints += XOR([a0, f[i]], Z[n-2-i])

        #9. (x_n-2-i_1, g_i) --AND--> (e_i)
        Constraints += AND([a1, g[i]], e[i])

        #10. (x_n-2-i_2, e_i) --XOR--> (v_i)
        Constraints += XOR([a2, e[i]], v[i])

    #11. (x_0, v_n-3) --XOR--> (z_0)
    Constraints += XOR([X[0], v[n-3]], Z[0])


    return (Variables, Constraints)


if __name__ == "__main__":
    n = 4
    addTrack = 1
    X = ['x0', 'x1', 'x2', 'x3']
    Y = ['y0', 'y1', 'y2', 'y3']
    Z = ['z0', 'z1', 'z2', 'z3']

    variables, constraints = ConstraintsByModAddVariables(X, Y, Z, 1)

    #variables, constraints = ConstraintsByModAddConstant(X, Z, 1)


    for c in constraints:
    	print (c)
