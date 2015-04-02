from sympy import Matrix
import sympy

def diag(x):
    def m(i, j):
        return x[i, i]
    try:
        a = Matrix(x.cols, 1, m)
    except:
        a = sympy.diag(*x)
    return a

def Trace(x):
    return sum(diag(x)

def W(H, A_inv):
    return (H.T*A_inv*H)**-1

def G(A_inv, H):
    return A_inv*H

def e(A_inv, f_F, H, beta_hat):
    return A_inv*(f_D - H*beta_hat)

