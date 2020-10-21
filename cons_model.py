import numpy as np
from scipy.linalg import polar
class Corotated:
    def __init__(self, mu, lambd, F):
        self.mu = mu
        self.lambd = lambd
        self.F = np.copy(F)
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
        ans[0,0] = 2*self.mu # +self.lambd*s1*s1
        ans[3,3] = 2*self.mu # +self.lambd*s0*s0
        ans[1,1] = ans[2,2] = 2*self.mu #- 2*self.mu/(s0+s1)
        ans[1,2] = ans[2,1] = -self.lambd*(self.J-1) #+ 2*self.mu/(s0+s1)
        ans[0,3] = ans[3,0] = self.lambd*(self.J - 1) # + self.lambd
        R11,R12,R21,R22 = self.R.flatten()
        dRdF = np.array([[R12*R12,-R12*R11,R12*R22,-R12*R21],
                         [-R11*R12,R11*R11,-R11*R22,R11*R21],
                         [R22*R12,-R22*R11,R22*R22,-R22*R21],
                         [-R21*R12,R21*R11,-R21*R22,R21*R21]])
        F11, F12, F21, F22 = self.F.flatten()
        dJdF = np.array([F22,-F21,-F12,F11])

        # return ans
        return ans-2*self.mu/(s0+s1)*dRdF+self.lambd*np.outer(dJdF,dJdF)

    def dP(self, dF):
        dR = polar(self.F+dF)[0]-self.R
        JFTinv = self.J*np.linalg.inv(self.F).transpose()
        dJFTinv = np.array([[dF[1,1], -dF[1,0]], [-dF[0,1], dF[0,0]]])
        return 2*self.mu*(dF-dR) + self.lambd*self.J*self.J*JFTinv*np.tensordot(JFTinv,dF)+ self.lambd*(self.J-1)*dJFTinv

if __name__ == "__main__":
    # Test P ~ dPdF
    mu, lambd = 1,1
    epsilon = []
    e = 0.1
    for i in range(5):
        epsilon.append(e)
        e /= 2
    print(epsilon)
    F = np.array([[1,1],[0,1]])
    print(np.linalg.inv(F).transpose())
    V = np.array([[1,0],[0,1]])
    model = Corotated(mu, lambd, F)
    DPDF = model.dPdF()
    print(DPDF)
    for i in range(5):
        print(f"Iteration {i}")
        model_plus = Corotated(mu,lambd,F+epsilon[i]*V)
        P_plus = model_plus.P()
        dP_plus = model.dP(epsilon[i]*V)
        model_minus = Corotated(mu,lambd,F-epsilon[i]*V)
        P_minus = model_minus.P()
        dP_minus = model.dP(-epsilon[i]*V)
        print(((P_plus-P_minus)/2/epsilon[i]).flatten())
        print(((dP_plus-dP_minus)/2/epsilon[i]).flatten())
        print(DPDF.dot(V.flatten()))
