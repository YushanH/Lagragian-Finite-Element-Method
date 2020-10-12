import numpy as np
from scipy.linalg import polar
class Corotated:
    def __init__(self, mu, lambd, F):
        self.mu = mu
        self.lambd = lambd
        self.F = F
        self.d = F.shape[0]
        self.U, self.Sigma, self.VT = np.linalg.svd(F)
        self.Sigma = np.diag(self.Sigma)
        self.R, self.S = self.U.dot(self.VT), self.VT.transpose().dot(self.Sigma).dot(self.VT)
        self.J = np.linalg.det(F)
        self.I = np.identity(self.d)

    def psi(self):
        ans = 0
        for i in range(self.d):
            ans += self.mu*(self.Sigma[i][i]-1)**2
        return ans + self.lambd/2*(self.J-1)**2

    def P(self):
        try:
            return 2*self.mu*(self.F - self.R) + self.lambd*(self.J-1)*self.J*np.linalg.inv(self.F).transpose()
        except:
            raise ValueError(self.F)
    def dPdF(self):
        ans = np.zeros((self.d*self.d, self.d*self.d))

        s0, s1 = self.Sigma[0,0], self.Sigma[1,1]
        ans[0,0] = 2*self.mu+self.lambd*s1*s1
        ans[3,3] = 2*self.mu+self.lambd*s0*s0
        ans[1,1] = ans[2,2] = 2*self.mu*(1-1/(s0+s1))
        ans[1,2] = ans[2,1] = 2*self.mu/(s0+s1)-self.lambd*(self.J-1)
        ans[0,3] = ans[3,0] = self.lambd*(2*self.J - 1)
        return ans



