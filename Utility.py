import numpy as np
import lmfit
from lmfit import  minimize, Parameters, report_fit
from scipy.optimize import curve_fit
from pathlib import Path

# for the numerical integration of mRNA and protein trafficking models
delta_x = 0.24
model_L = 500
# possbile lengths for dendritic distributions
Lengths = np.asarray([2,25,50,75,100,125,150,200,250])
bin_size = 5
bins = np.arange(0, Lengths.max(), bin_size)
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
COLORS_dict = {"spine":"#005f73","shaft":'#CA6702',"spine_s":"#005f73","spine_i":'#CA6702',"shaft_s":"#005f73","shaft_i":'#CA6702'}
def Area(radius):
    return np.pi*(radius**2)

def CreateFolderRecursive(folder):
    Path(folder).mkdir(parents=True, exist_ok=True)
def GetRNAPunctaStat(puncta):
    # the puncta is of the form (x,y,r,max,min,mean.std,median, insoma,distance_from_soma)
    x, y, r, mx, mn, mu, delta, med, inSoma, dfs = puncta

    area = Area(r)
    stat1 = area * mu
    stat2 = area * med
    stat3 = mu
    stat4 = med
    stat5 = area
    return np.array([stat1, stat2, stat3, stat4, stat5, 1, dfs])
def GetDictSum(d):
    dict_sum = np.zeros(6)
    for key in d.keys():
        dict_dum = np.sum(dict_sum, GetRNAPunctaStat(d[key]))

def GetMatSum(mat,ax=0):
    return np.sum(mat,axis=ax)

def SortPunctas(ps,column=0):
    return ps[ps[:,column].argsort()]


def ChiSq(yd,y_fit,sigmas):
    nzs = np.nonzero(sigmas)
    # print(nzs)
    r_yd = np.take(yd,nzs)
    r_yf = np.take(y_fit,nzs)
    r_sgs = np.take(sigmas,nzs)
    # print("r_sgs = ",r_sgs)
    residuals = r_yd - r_yf
    chi_squ = np.sum((residuals/ r_sgs)**2)
    # print(residuals,r_sgs)
    return chi_squ

def BinnedSum(arr, bins, num_col=-1, name=None):
    # print("binned_sum for =",name)

    if len(arr.shape) == 2:
        rr, cc = arr.shape
        binned_sum = np.zeros((len(bins), cc))
        arr = SortPunctas(arr, num_col)
        # max_lenght = arr[-1,-1]
        digitized = bins.searchsorted(arr[:, num_col])

        #     breakpoint()
        for c in range(0, cc):
            try:
                binned_sum[:, c] = np.bincount(digitized, weights=arr[:, c], minlength=len(bins))
            except:
                print("puncta is ", arr)
                # breakpoint()
        binned_sum[:, num_col] = bins
        return binned_sum
    else:
        print("quite not the shape", arr.shape)
        return np.zeros((len(bins), 6))


def GetUniqueRows(mat):
    return np.unique(mat, axis=0)

def ExpFit2(xdata,ydata,sigmas,Fidx,Lidx,molecule):
    pars = Parameters()
    pars.add('amplitude',1,vary=False)
    pars.add('decay',1,min=0)
    mod = lmfit.models.ExponentialModel()
    out = mod.fit(ydata[Fidx:], pars, x=xdata[Fidx:])
    y_fit = exponential(xdata[Fidx:],-1.0/out.params['decay'])
    residuals = ydata[Fidx:]- y_fit[Fidx:]
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata[Fidx:]-np.mean(ydata[Fidx:]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    chi_squ = np.sum((residuals/sigmas)**2)
    print("here",chi_squ)
    # breakpoint()
    return y_fit,r_squared,out.chisqr


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def normExponential(x, params):
    b = params['b'].value
    return np.exp(b * x)


def oneExponential(x, params):
    a = params['a'].value
    b = params['b'].value
    return a * np.exp(b * x)


def twoExponential(x, params):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    return a * np.exp(b * x) + c * np.exp(d * x)


def ExpFit(ftype, xdata, ydata, sigmas, Fidx, Lidx, molecule):
    """
        Fit a function to a given distribution

    """
    if ftype == "NormE":
        param_bounds = ([-np.inf], [np.inf])
        popt, pcov = curve_fit(normExponential, xdata[Fidx:], ydata[Fidx:], bounds=param_bounds, maxfev=5000)
        y_fit = normExponential(xdata, *popt)
    elif ftype == "1E":
        param_bounds = ([-np.inf, -np.inf], [+np.inf, 0])
        popt, pcov = curve_fit(oneExponential, xdata[Fidx:], ydata[Fidx:], bounds=param_bounds)
        y_fit = oneExponential(xdata, *popt)
    elif ftype == "2E":
        param_bounds = ([-np.inf, -np.inf, -np.inf, -np.inf], [+np.inf, 0, +np.inf, 0])
        popt, pcov = curve_fit(twoExponential, xdata[Fidx:], ydata[Fidx:], bounds=param_bounds)
        y_fit = twoExponential(xdata, *popt)
    else:
        raise NotImplementedError("ftype: {} not implemented, contact author or define it yourself".format(ftype))

    print("fitted " + ftype, popt)
    residuals = ydata[Fidx:] - y_fit[Fidx:]
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata[Fidx:] - np.mean(ydata[Fidx:])) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    chi_squ = ChiSq(ydata[Fidx:], y_fit[Fidx:], sigmas[Fidx:])
    print("chi-squared = ", chi_squ)
    return y_fit, r_squared, chi_squ


def ExpFitWithMinimize(ftype, xdata, ydata, sigmas, Fidx, Lidx, molecule):
    """
        Fit a function to a given distribution using lmfit minimize method

    """
    fit_paramas = Parameters()
    # np.random.seed(2022)
    exp_min = -2.0
    exp_max = 0
    pref_min = 0
    pref_max = 200

    if ftype == "NormE":
        b_init = np.random.uniform(exp_min, exp_max)
        fit_paramas.add('b', b_init, min=exp_min, max=exp_max)

        out2 = minimize(Residual, params=fit_paramas, method='leastsq',
                        args=(normExponential, xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = normExponential(xdata[Fidx:], out2.params)
        chi_squ = ChiSq(ydata[Fidx:], y_fit, sigmas) / (ydata[Fidx:].shape[0] - 1)
        # breakpoint()

    elif ftype == "1E":

        a_init = np.random.uniform(pref_min, pref_max)
        b_init = np.random.uniform(exp_min, exp_max)
        fit_paramas.add('a', a_init, min=pref_min, max=pref_max)
        fit_paramas.add('b', b_init, min=exp_min, max=exp_max)
        out2 = minimize(Residual, params=fit_paramas, method='leastsq',
                        args=(oneExponential, xdata[Fidx:], ydata[Fidx:]))
        print("reporting 1E fits ")
        report_fit(out2.params)
        y_fit = oneExponential(xdata[Fidx:], out2.params)
        chi_squ = ChiSq(ydata[Fidx:], y_fit, sigmas) / (ydata[Fidx:].shape[0] - 2)
    elif ftype == "2E":
        print("fitting 2E")
        a_init = np.random.uniform(pref_min, pref_max)
        b_init = np.random.uniform(exp_min, exp_max)
        c_init = np.random.uniform(pref_min, pref_max)
        d_init = np.random.uniform(exp_min, exp_max)
        fit_paramas.add('a', a_init, min=pref_min, max=pref_max)
        fit_paramas.add('b', b_init, min=exp_min, max=exp_max)
        fit_paramas.add('c', c_init, min=pref_min, max=pref_max)
        fit_paramas.add('d', d_init, min=exp_min, max=exp_max)
        out2 = minimize(Residual, params=fit_paramas, method='leastsq',
                        args=(twoExponential, xdata[Fidx:], ydata[Fidx:]))
        print("reporting 2E fits ")
        report_fit(out2.params)
        y_fit = twoExponential(xdata[Fidx:], out2.params)
        chi_squ = ChiSq(ydata[Fidx:], y_fit, sigmas) / (ydata[Fidx:].shape[0]-4)
    else:
        raise NotImplementedError("ftype: {} not implemented, contact author or define it yourself".format(ftype))

    # print("fitted "+ ftype, popt)



    return y_fit, chi_squ


def Residual(paras, fun, x, data):
    expected_vals = fun(x, paras)
    res = expected_vals - data
    return res

def getminmax(arr,orig_min=np.Inf,orig_max=-np.Inf):
    if np.min(arr) < orig_min:
        orig_min = np.min(arr)
    if np.max(arr) > orig_max:
        orig_max = np.max(arr)
    return orig_min,orig_max


def exp_fit(x,a,b):
    return a*np.exp(-b*x)

def line_fit(x, m,c):
    return m*x + c
def R_seq(y_fit,y_orig):
    ss_res = ((y_orig-y_fit)**2).sum()
    ss_tot = ((y_orig-y_orig.mean())**2).sum()
    # print("in R_Seq =",ss_tot,ss_res)
    return 1 - (ss_res/ss_tot)