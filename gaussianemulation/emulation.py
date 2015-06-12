from sympy import *
from mpmath import *
from util import *

def norm_points(points, high=None, low=None):
    if high is None:
        high = max(points)
    if low is None:
        low = min(points)
    return Matrix((points-low)/(high - low))

def cor(x, x_prime, delta):
    """
    Calculate the correlation between two points in n-dimensional
    space.
    """
    def m(i, j):
        #print x[i, j], x_prime[i, j], delta[i]
        return ((x[i, j]-x_prime[i, j])/delta[i])**2
    return exp(-sum(Matrix(x.rows, x.cols, m)))

def cov(x, x_prime, delta, sigma_squared):
    return sigma_squared*cor(x, x_prime, delta)

def pi_star_delta(sigma_hat_sq, n, q, A_inv, A, H):
    return power(sigma_hat_sq, -(n-q)/2.)*power(A.det(), -1/2.)*power((H.T*A_inv*H).det(), -1/2.)

def sigma_hat_sq(n, q, f_D, A_inv, H):
    return (((n-q-2)**-1.)*f_D.T*(A_inv-A_inv*H*((H.T*A_inv*H)**-1)*H.T*A_inv)*f_D)[0, 0]

def sigma_hat_sq_beta(n, q, f_D, A_inv, H_, beta_hat):
    return ((f_D-H_*beta_hat).T*A_inv*(f_D-H_* beta_hat))/(n-q-2.)

def A(points, delta):
    """
    Autocorrelation matrix.
    """
    def m(i, j):
        #print "next", points[:,i], points[:,j], delta
        return cor(points[:,i], points[:,j], delta)
    a = Matrix(points.cols, points.cols, m)
    return a

def fast_A(points, delta):
    scales = diag(map(lambda x: mpf(1)/x**2, delta))
    #print points, scales
    R = points.T * scales * points
    S = Matrix([diag(R).T] * R.cols)
    A = (R + R.T - S - S.T).applyfunc(exp)
    return A

def tau(delta):
    return 2*log(delta)

def pi_star_tau(n, q, f_D, A_inv, A, H):
    return pi_star_delta(sigma_hat_sq(n, q, f_D, A_inv, H), n, q, A_inv, A, H)

def H(x, h):
    return h(x).T

def beta_hat(H, A_inv, f_D):
    return ((H.T*A_inv*H)**-1)*H.T*A_inv*f_D

def h(x):
    return Matrix(ones(1, x.cols)).T.col_join(x)

def mahalanobis(f_D, m_star, v_star):
    """
    Calculate the mahalanobis distance for a trained emulator.
    """
    return (f_D - m_star).T*(v_star**-1)*(f_D - m_star)

def se(f_D, m_star, v_star):
    return (f_D - m_star).T / diag(v_star).applyfunc(sqrt)

def make_t(D, delta_hat):
    """
    Make a t (c, on the MUCM toolkit) function.
    """
    def f(x):
        def m(i, j):
            return cor(D[:,i], x[:,j], delta_hat)
        a = Matrix(D.cols, x.cols, m)
        return a
    return f

def make_m_star(h, b_hat, A_inv, f_D, t):
    def f(x):
        return h(x).T*b_hat + t(x).T*A_inv*(f_D - H_*b_hat)
    return f

def make_v_star(d_hat_sq, delta_hat, A_inv, h, H_, t):    
    """
    Make the posterior variance function.
    """
    def f(x, x_prime):
        def m(i, j):
            return cor(x[:,i], x_prime[:,j], delta_hat)
        c_x_x_prime = Matrix(x.cols, x_prime.cols, m)
        return d_hat_sq*(c_x_x_prime - t(x).T*A_inv*t(x_prime) + (h(x).T-t(x).T*A_inv*H_)*((H_.T*A_inv*H_)**-1)* (h(x_prime).T-t(x_prime).T*A_inv*H_).T)
    return f