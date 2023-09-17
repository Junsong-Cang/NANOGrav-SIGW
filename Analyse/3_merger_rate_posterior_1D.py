from src.merger import *

reload = 1
FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
ResultFile='/Users/cangtao/Desktop/PBH_posteriors.pdf'
LineWidth = 2
FontSize = 18

import matplotlib.pyplot as plt
import getdist, os
from getdist import plots

RateRoot = FileRoot+'Merger_Rate'

RateTab = np.loadtxt(FileRoot +'Merger_Rate.txt')
lf_vec = RateTab[:,0]
lm_vec = RateTab[:,1]
sbh_vec = RateTab[:,0]

def model(theta):
    lf, lm, sbh = theta[5], theta[6], theta[7]
    dist = (lf_vec - lf)**2 + (lm_vec - lm)**2 + (sbh_vec - sbh)**2
    idx = np.argmin(dist)
    r = np.log10(RateTab[idx, 3])
    return r

if reload:
    t1 = PyLab.TimeNow()
    PyLab.derived_param_chains(
        model = model,
        old_root = FileRoot,
        new_root = FileRoot + 'Merger_Rate_1D',
        derived_names = [{'name': 'R', 'latex': '\log_{10}R_{z0}'}],
        clean_up = 0,
        show_status = 1)
    PyLab.Timer(t1)
