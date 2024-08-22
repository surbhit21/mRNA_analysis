import json
import matplotlib
import numpy as np

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from Utility import *
class Protein_model():
    """
        Class for the steady state solution of protein model equations from Fonkeu et al. 2019 under closed BC at x = L
        and constant flux at x=0:

    """

    def __init__(self, D_P, v_P, half_life_protein, JPin,dx,L):
        self.D_P = D_P
        self.v_P = v_P
        self.t_half_protein = half_life_protein
        self.k_P = np.log(2) / (self.t_half_protein * 24 * 60 * 60);
        self.JPin = JPin
        self.dx = dx
        self.L = L
        self.x_grid = np.arange(0, self.L, dx)

    # def sanity_check(self):
    #     assert self.r_dist.shape == self.x_grid.shape       # the mRNA distribution array should be defined on each grid point

    def GetSSDist(self):
        lambda_r1 = (-self.v_P + np.sqrt(self.v_P ** 2 + (4 * self.D_P * self.k_P))) / (2 * self.D_P)
        lambda_r2 = (-self.v_P - np.sqrt(self.v_P ** 2 + (4 * self.D_P * self.k_P))) / (2 * self.D_P)
        C1 = (self.JPin * np.exp(-1. * lambda_r2 * self.L)) / (
                    (self.v_P * self.D_P * lambda_r1) * (np.exp(-1. * lambda_r1 * self.L) - np.exp(-1. * lambda_r2 * self.L)))
        C2 = (self.JPin * np.exp(-1. * lambda_r1 * self.L)) / (
                    (self.v_P * self.D_P * lambda_r2) * (np.exp(-1. * lambda_r2 * self.L) - np.exp(-1. * lambda_r1 * self.L)))
        print(lambda_r1, lambda_r2, C1, C2)

    def fun_protein(self, x, y):
        """
        function that defines the darivatives
        """
        p, dp = y
        # print(y.shape)
        return [
            dp,
            (self.v_P * dp + self.k_P * p) / self.D_P
        ]

    def bc_protein(self, ya, yb):
        """
            function for boundary condition values
        """
        return np.array(
            [ya[1] - self.v_P * ya[0] / self.D_P + self.JPin / self.D_P, yb[1] - self.v_P * yb[0] / self.D_P])

    def SolveNumericalProtein(self):
        y = np.zeros((2, self.x_grid.size))
        soln = solve_bvp(self.fun_protein, self.bc_protein, self.x_grid, y, max_nodes=int(self.L/self.dx), verbose=0, tol=1e-3, bc_tol=1e-8)
        # breakpoint()
        p_dist = soln.y[0]
        dp_dist = soln.y[1]
        self.x_grid = soln.x
        # breakpoint()
        # print(self.bc(soln.y[:,0],soln.y[:,1]))
        return self.x_grid, p_dist




def RunSSProtein(D_P=0.22, v_P=1.4e-3,t_half = 3.12,x_range = [0,500]):
    Jpin = 0.01
    dx = 1
    # r_dist = oneExponential(x_grid,fit_paramas)
    protein_ss = Protein_model(D_P, v_P, t_half, Jpin,dx,500)
    # protein_ss.sanity_check()
    breakpoint()
    x_grid, p_dist = protein_ss.SolveNumericalProtein()
    param_dict = {}
    param_dict["D_P"] = protein_ss.D_P
    param_dict["v_P"] = protein_ss.v_P
    param_dict["half_life_protein"] = protein_ss.t_half_protein
    param_dict["dx"] = protein_ss.dx
    param_dict["JPin"] = protein_ss.JPin
    # breakpoint()
    # with open("./TotalAMPAModelParams.json", "w") as fp:
    #     json.dump(param_dict, fp)
    # fp.close()
    # plt.plot(protein_ss.x_grid,p_dist,label="Glua2")
    # plt.legend()
    # plt.show()
    return x_grid[x_range[0]//dx:x_range[-1]//dx], p_dist[x_range[0]//dx:x_range[-1]//dx]

# RunSSProtein(0.61,1e-5)
def RunModelWithFile(param_file):
    with open(param_file, "r") as fp:
        params = json.load(fp)
    fp.close()
    CNIH2_model1 = Protein_model(**params);
    CNIH2_model1.r_dist= np.asarray(CNIH2_model1.r_dist)
    # CNIH2_model1.x_grid = np.asarray(CNIH2_model1.x_grid)
    CNIH2_model1.sanity_check()
    x_grid,CNIH2_dist = CNIH2_model1.SolveNumericalProtein()
    return CNIH2_dist, CNIH2_model1

# dist,model_obj = RunModelWithFile("./CNIH2SSModelParamsTemporal.json")
# breakpoint()