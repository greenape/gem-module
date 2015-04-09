from sympy import *
from mpmath import *
from util import *
from emulation import *


def fit(n, q, f_D, D, h, start):
    def f(x):
        x = Matrix(map(lambda y: exp(mpf(y)/2.), x))
        a = A(D, x)
        a_inv = a**-1
        return -log(pi_star_tau(n, q, f_D, a_inv, a, H(D, h)))
    tau = minimize(f, start, method="nelder-mead", options={'disp':True})['x']
    delta_hat = map(exp, map(lambda x: x/2., map(mpf,tau)))
    return delta_hat

class Emulator(object):
    """
    A gaussian emulator.
    """
    def __init__(self, inputs, outputs):
        self.inputs = Matrix(inputs).T
        self.outputs = Matrix(outputs)
        self.n = self.inputs.cols
        self.p = self.inputs.rows
        self.q = self.p + 1

    def fit(self):
        """
        Attempt to fit the emulator, given the inputs and outputs
        provided.
        """
        # Fitted emulator
        self.A_inv = A(D, self.delta_hat)**-1
        self.beta_hat = beta_hat(H(self.points, h), self.A_inv, self.outputs)
        self.sigma_hat_sq = sigma_hat_sq(self.n, self.q, self.outputs, self.A_inv, H(self.points, h))
        self.H = H(self.points, h)
        self.t = make_t(self.points, self.delta_hat)
        self.m_star = make_m_star(h, self.beta_hat, self.A_inv, self.outputs, self.t)
        self.v_star = make_v_star(self.sigma_hat_sq, self.delta_hat, self.A_inv, h, self.H, self.t)

    def validate(self, D_prime, f_D_prime):
        m_star_D_prime, v_star_D_prime = self.predict(D_prime)
        self.ses = se(f_D_prime, m_star_D_prime, v_star_D_prime)
        self.mahalanobis = mahalanobis(f_D_prime, m_star_D_prime, v_star_D_prime)

    def predict(self, D_prime):
        """
        Predict new points using the fitted emulator.
        """
        try:
            m_star_D_prime = self.m_star(D_prime)
            v_star_D_prime = self.v_star(D_prime, D_prime)
            return m_star_D_prime, v_star_D_prime
        except AttributeError:
            raise

    def add_points(self, D_prime, f_D_prime):
        """
        Add additional points to the emulator, and refit.
        """
        self.inputs.row_join(D_prime)
        self.outputs.col_join(f_D_prime)
        self.fit()
