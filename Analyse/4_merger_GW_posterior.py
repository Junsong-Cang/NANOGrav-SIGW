from src.merger import *

reload = 0
FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'
v = np.logspace(-3, 10, 80)
LineWidth = 2
FontSize = 16
ncpu = 12
nm = 50
nz = 50

import matplotlib.pyplot as plt
import os

LogFile_0 = FileRoot + 'GW_0.txt'
LogFile_1 = FileRoot + 'GW_1.txt'

def model_0(theta):
    Use_S2 = 0
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    g = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh, 
        v = v,
        Fast = 0,
        mf_model = 0, 
        sbh_width = 6, 
        nm = nm,
        nz = nz,
        show_status = 0,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    PyLab.SaySomething('tmp_2.txt')

    # Print results to a file
    
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    nv = len(v)
    for idx in np.arange(0, nv):
        Str = Str + '  {0:.5E}'.format(g[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()

    return g

def model_1(theta):
    Use_S2 = 1
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    g = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh, 
        v = v,
        Fast = 0,
        mf_model = 0, 
        sbh_width = 6, 
        nm = nm,
        nz = nz,
        show_status = 0,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    PyLab.SaySomething('tmp_2.txt')
    
    # Print results to a file
    
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    nv = len(v)
    for idx in np.arange(0, nv):
        Str = Str + '  {0:.5E}'.format(g[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()
    
    return g

if reload:
    try:
        os.remove(LogFile_0)
    except:
        pass
    try:
        os.remove(LogFile_1)
    except:
        pass
    t1 = PyLab.TimeNow()
    r0 = PyLab.mcmc_derived_stat(model_function = model_0, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    r1 = PyLab.mcmc_derived_stat(model_function = model_1, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    # np.savez('data/3_merger_GW_posterior.npz', r0 = r0, r1 = r1)
    PyLab.Timer(t1)
