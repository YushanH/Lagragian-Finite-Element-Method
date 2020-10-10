import numpy as np
import scipy
class Corotated_element:
    def __init__(self, mu, lam, F):
        self.mu = mu
        self.lambd = lambd
        self.F = F
        self.d = F.shape[0]
        self.R, self.S = scipy.linalg.polar(F)
        self.I = np.identity(self.d)
    def self.psi():
        norm = np.linalg.norm(self.S-self.I)
        trace = np.trace(self.S-self.I)
        return self.mu*norm**2 + self.lambd/2*trace**2

    def self.P():
        return 2*self.mu*(self.F - self.R) + self.lambd*np.trace(np.transpose(self.R)*self.F-self.I)*self.R

    def self.dPdF():
        pass


