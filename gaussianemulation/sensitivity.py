from sympy import *
from mpmath import *
from util import *

def posterior_mean_i(r_i, r_h, t_i, r_t, b_hat, e_):
    return (r_i - r_h.T)*b_hat + (t_i - r_t.T)*e_

def posterior_mean_v_w(sigma, beta, A_inv, e_, w, H, u, u_w, p_w, q_w, r, t):
    """
    Sensitivity variance
    """
    term = q_w - s_w*A_inv*H
    term -= H.T*A_inv*s_w.T
    term += H.T*A_inv*p_w*A_inv*H
    
    result = u_w-Trace(A_inv*p_w) 
    result += Trace(w*term)
    result *= sigma
    result += (e_.T*p_w*e_)[0]
    result += (2*beta.T*s_w*e_)[0]
    result += (beta.T*q_w*beta)[0]
    #second term
    rtah = r-t*A_inv*H
    result -= (sigma*(u-t*A_inv*t.T + rtah*w*rtah.T))[0]
    result -= power((r*beta + t*e_)[0], 2.)
    return result

def R_w(m, x, w_bar, w):
    """
    1xq vector
    """
    #return np.hstack(([1], m.T[0, w_bar], x.T[0, w]))
    #return np.hstack(([1], x.T[0, w], m.T[0, w_bar]))
    indices = w+w_bar
    indices.sort()
    return Matrix([1]).row_join(Matrix(map(lambda k: x.T[0, k] if k in w else m.T[0, k], indices)).T)

def T_w(xw, x, b, c, m, w, w_bar, v):
    return Matrix(1, x.cols, lambda foo, k: T_w_element(xw, x, B, C, m, w, w_bar, v, k))

def T_w_element(xw, x, b, c, m, w, w_bar, v, k):
    """
    1xn vector, R_t is a special case for w is empty
    """
    #print xw, x, b, c, m, w, w_bar, v, k
    exp2 = exp(-0.5*(xw[w]-x[w,k])*2*c[w,w]*(xw[w]-x[w,k]))
    def exp1(i):
        return exp(-0.5*((2*c[i, i]*b[i, i])/(2*c[i, i] + b[i, i]))*power(x[i,k]-m[i], 2.))
    def div(i):
        return sqrt(b[i,i])/sqrt(2*c[i,i] + b[i,i])
    return (1 - v)*reduce(mul, map(lambda i: div(i)*exp1(i)*exp2, w_bar))


def posterior_mean_m_w(r_w, b_hat, t_w, e_):
    return r_w*b_hat + t_w*e_

## Sensitivity special case 2

def U_w(v, b, c, w_bar):
    return (1 - v)*reduce(mul, map(lambda i: power(b[i, i]/(b[i, i]+2*(2*c[i, i])), 0.5), w_bar))

def P_w(x, m, b, c, w):
    return Matrix(x.cols, x.cols, lambda k, l: P_w_element(x, k, l, m, B, C, w))

def P_w_element(x, k, l, m, b, c, w):
    """
    nxn matrix
    """
    def term_1(i):
        return (b[i, i]/(2*c[i, i] + b[i, i]))*exp(-0.5*((2*c[i, i]*b[i, i])/(2*c[i, i]+b[i, i]))*(power(x[i, k] - m[i], 2.) + power(x[i, l]-m[i], 2.)))
    def term_2(i):
        return power(b[i,i]/(4*c[i,i] + b[i, i]), 0.5)*exp(-0.5*(1/(4*c[i, i] + b[i, i]))*(4*power(c[i,i], 2.)*power(x[i,k]-x[i,l],2.0) + 2*c[i,i]*b[i,i]*(power(x[i,k]-m[i], 2.) + power(x[i,l]-m[i], 2.))))
    return power(1 - v, 2.)*reduce(mul, map(term_1, w_bar))*reduce(mul, map(term_2, w))

def P(r_t):
    return r_t.T*r_t

def S_w(w, w_bar, b, c, x, m, v):
    return Matrix(x.cols, q, lambda l, k: S_w_element(k, l, w, w_bar, B, C, x, m, v)).T

def S_w_element(k, l, w, w_bar, b, c, x, m, v):
    """
    qxn matrix
    """
    h_k = 1
    k -= 1
    if k in w_bar:
        h_k = m[k]
    elif k in w:
        h_k = (2*c[k,k]*x[k, l] + b[k, k]*m[k])/(2*c[k, k] + b[k, k])
    def t2(i):
        return power(b[i, i], .5) / power(2*c[i,i] + b[i, i], 0.5)*exp(-0.5*(((2*c[i, i]*b[i,i])/(2*c[i, i] + b[i, i]))*power(x[i, l] - m[i], 2.)))
    return (1 - v)*h_k*reduce(lambda x, y: x*y, map(t2, w + w_bar))

def S(r_h, r_t):
    return r_h.T*r_t

def Q_w(m, w, w_bar, b):
    """
    qxq matrix
    """
    indices = w+w_bar
    indices.sort()
    
    result = Matrix([1]).row_join(m[w_bar,0].T).row_join(m[w,0].T)
    result = result.col_join(m[w_bar,0].row_join(m[w_bar,0]*m[w_bar,0].T).row_join(m[w_bar,0]* m[w,0].T))
    result = result.col_join(m[w,0].row_join(m[w,0]* m[w_bar,0].T).row_join(m[w,0]*m[w,0].T + (B**-1)[w, w]))
    if w[0] < w_bar[0]:
        x = result.col_insert(1,result[:,result.cols-1])
        x.col_del(x.cols-1)
        result = x.row_insert(1, x[x.rows-1,:])
        result.row_del(result.rows-1)
    return result
    
def Q(r_h):
    return r_h.T*r_h