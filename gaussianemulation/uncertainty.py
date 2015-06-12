from sympy import *
from mpmath import *
from util import *

def E(r_h, b_hat, r_t, e):
    return r_h.T*b_hat + r_t.T*e

def V(sigma, u, r_t, A_inv, r_h, g, w):
    res = Matrix([u])
    res -= r_t.T*A_inv*r_t
    res += (r_h - g.T*r_t).T*w*(r_h-g.T*r_t)
    res *= sigma
    return sigma*(u - r_t.T*A_inv*r_t + (r_h - g.T*r_t).T*w*(r_h-g.T*r_t))

def do_E_var(i_i, i_2, V_, E_):
    return (i_i-v) + (i_2 - power(E_, 2.))

def E_var():
    r_tt = R_tt(D, C, B, m, v)
    r_hh = R_hh(m, B)
    r_ht = R_ht(D, B, C, v, m, h)
    i_1 = I_1(s_hat_sq, A_inv, r_tt, w, r_hh, r_ht, g)
    i_2 = I_2(b_hat, r_hh, r_ht, e_, r_tt)
    return do_E_var(i_1, i_2[0,0], V_[0,0], E_[0,0])

def I_1(sigma, A_inv, r_tt, w, r_hh, r_ht, g):
    return sigma*(mpf(1)-Trace(A_inv*r_tt) + Trace(w*(r_hh - 2*r_ht*g + g.T*r_tt*g)))

def I_2(beta, r_hh, r_ht, e_, r_tt):
    return beta.T*r_hh*beta + 2*beta.T*r_ht*e_ + e_.T*r_tt*e_

def Q_kl(x, xk, xl, C, B, m):
    return 2*(x - xk).T*C*(x - xk) + 2*(x - xl).T*C*(x - xl) + (x - m).T*B*(x - m)

def Q_k(x, xk, m, B, C):
    return (2*(x - xk).T*C*(x - xk) + (x-m).T*B*(x-m))[0,0]

def m_kl(xk, xl, C, B, m):
    return ((4*C + B)**-1)*(2*C*xk + 2*C*xl + B*m)

def m_k(x, C, B, m):
    return ((2*C + B)**-1)*(2*C*x + B*m)

def R_h(m):
    return Matrix([1]).col_join(m)

def R_hh(m, B):
    #np.vstack((np.hstack(([[1]], m.T)), np.hstack((m, m.dot(m.T) + B.getI()))))
    return Matrix([1]).row_join(m.T).col_join(m.row_join(m*m.T + B**-1))

def R_ht(D, B, C, v, m, h):
    return reduce(lambda x, y: x.row_join(y),map(lambda k: R_ht_elem(D, k, B, C, v, m, h), range(D.cols))) #matrix

def R_ht_elem(X, k, B, C, v, m, h):
    x = X[:,k]
    m_prime_k = m_k(x, C, B, m)
    return R_t(X, k, B, C, v, m)*Matrix([1]).col_join(m_prime_k)

def R_tt(D, C, B, m, v):
    return Matrix(D.cols, D.cols, lambda i, j: R_tt_element(D, i, j, C, B, m, v))

def R_tt_element(x, k, l, C, B, m, v):
    xk = x[:,k]
    xl = x[:,l]
    qkl = Q_kl(m_kl(xk, xl, C, B, m), xk, xl, C, B, m)[0,0]
    return power(1-v, 2.)*power(det(B), 0.5)*power(det(4*C + B), -0.5)*exp(- qkl/2.)

def R_t(D, B, C, v, m):
    return Matrix(map(lambda k: R_t_elem(D, k, B, C, v, m), range(D.cols)))

def R_t_elem(X, k, B, C, v, m):
    X = X[:,k]
    m_prime_k = m_k(X, C, B, m)
    q_k = Q_k(m_prime_k, X, m, B, C)
    
    return (1-v)*power(det(B), 0.5)*power(det(2*C + B), -0.5)*exp(-q_k/2.)