import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from mRNA_SS_model import *
from scipy.integrate import solve_bvp
from Utility import *
class Protein_model():
    """
        Class for the steady state solution of protein model equations from Fonkeu et al. 2019 under closed BC at x = L
        and constant flux at x=0:

    """

    def __init__(self, D_P, v_P, half_life_protein, JPin, beta_p,r_dist,dx):
        self.D_P = D_P
        self.v_P = v_P
        self.t_half_mRNA = half_life_protein
        self.k_P = np.log(2) / (self.t_half_mRNA * 24 * 60 * 60);
        self.JPin = JPin
        self.r_dist = r_dist
        self.beta_p = beta_p
        self.dx = dx
        self.x_grid = np.arange(0, model_L, dx)
        assert self.r_dist.shape == self.x_grid.shape       # the mRNA distribution array should be defined on each grid point

    def GetSSDist(self):
        lambda_r1 = (-self.v_P + np.sqrt(self.v_P ** 2 + (4 * self.D_P * self.k_P))) / (2 * self.D_P)
        lambda_r2 = (-self.v_P - np.sqrt(self.v_P ** 2 + (4 * self.D_P * self.k_P))) / (2 * self.D_P)
        C1 = (self.JPin * np.exp(-1. * lambda_r2 * model_L)) / (
                    (self.v_P * self.D_P * lambda_r1) * (np.exp(-1. * lambda_r1 * model_L) - np.exp(-1. * lambda_r2 * model_L)))
        C2 = (self.JPin * np.exp(-1. * lambda_r1 * model_L)) / (
                    (self.v_P * self.D_P * lambda_r2) * (np.exp(-1. * lambda_r2 * model_L) - np.exp(-1. * lambda_r1 * model_L)))
        print(lambda_r1, lambda_r2, C1, C2)

    def fun_protein(self, x, y):
        """
        function that defines the darivatives
        """
        p, dp = y
        # print(y.shape)
        return [
            dp,
            (self.v_P * dp + self.k_P * p - self.beta_p * self.r_dist[0:dp.shape[0]]) / self.D_P
        ]

    def bc_protein(self, ya, yb):
        """
            function for boundary condition values
        """
        return np.array(
            [ya[1] - self.v_P * ya[0] / self.D_P + self.JPin / self.D_P, yb[1] - self.v_P * yb[0] / self.D_P])

    def SolveNumericalProtein(self):
        y = np.zeros((2, self.x_grid.size))
        soln = solve_bvp(self.fun_protein, self.bc_protein, self.x_grid, y, max_nodes=int(model_L/delta_x), verbose=0, tol=1e-3, bc_tol=1e-8)
        # breakpoint()
        p_dist = soln.y[0]
        dp_dist = soln.y[1]
        self.x_grid = soln.x
        # breakpoint()
        # print(self.bc(soln.y[:,0],soln.y[:,1]))
        return self.x_grid, p_dist




def RunSSProtein(D_P, v_P):
    x_grid,r_dist = RunSS(0.23,0)
    protein_ss = Protein_model(D_P, v_P, 2, 0.001, 0.001,r_dist,delta_x)
    x_grid, p_dist = protein_ss.SolveNumericalProtein()
    return x_grid, p_dist
    # plt.plot(protein_ss.x_grid, r_dist / r_dist[0],label="mRNA")
    # plt.plot(protein_ss.x_grid,p_dist/p_dist[0],label="Protein")
    # plt.legend()
    # plt.show()
#
# RunSSProtein(10**0.9999996,10**-2.45334484)