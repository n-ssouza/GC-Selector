import numpy as np
from numpy.linalg import inv 

import scipy.sparse.linalg as sp


def MLPCA(X, Xsd, p, MaxIter=1e5):
    epsilon = 1e-10
    MaxIter = MaxIter
    m = X.shape[0]
    n = X.shape[1]
    VarX = np.multiply(Xsd, Xsd)
    U, o, V = sp.svds(X, k=p)
    i = 0
    Sold = 0
    k = -1

    while (k < 0):
        i = i + 1
        print('Iteration',i)
        Sobj = 0
        LX = np.matrix(np.zeros((m,n)))

        for j in range (0, n):
            Q = np.diagflat(1 / VarX[:, j])
            #F = inv(U.T @ Q @ U)
            F = np.linalg.inv(np.matmul(np.matmul(U.T, Q), U))
            #LX[:, j] = U @ (F @ (U.T @ (Q @ X[:, j])))
            LX[:, j] = np.matmul(U, (np.matmul(F, (np.matmul(U.T, (np.matmul(Q, X[:, j]))))))).reshape((m,1))
            Dx = np.matrix(X[:, j] - LX[:, j])
            #Sobj = Sobj + Dx.T @ Q @ Dx
            Sobj = Sobj + np.matmul(np.matmul(Dx.T, Q), Dx)
            print('Dim.',j,'DONE')
        
        if i % 2 == 1:
            ConvCalc = np.abs(Sold - Sobj) / Sobj
            if (ConvCalc < epsilon).any():
                k = 0
            if i > MaxIter:
                k = 1
                exit("MaxIter exceeded")

        if k < 0:
            Sold = Sobj
            U, o, V = sp.svds(LX, k=p)
            V = V.T
            X = X.T
            VarX = VarX.T
            n = X.shape[1]
            U = V

    U, o, V = sp.svds(LX, k=p)
    S = np.matrix(np.diag(o))

    return U, S, V

