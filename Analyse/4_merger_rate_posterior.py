from src.merger import *

reload = 1

FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'

LineWidth = 2
FontSize = 18
nm = 100
S_Tab_Len = 200
ncpu = 12
z = np.logspace(0, 2.2, 40) - 1

LogFile_0 = FileRoot + 'Merger_Rate_0.txt'
LogFile_1 = FileRoot + 'Merger_Rate_1.txt'

import matplotlib.pyplot as plt
import os

nz = len(z)

def model_0(theta):
    Use_S2 = 0
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    r = Merger_Rate(
        fbh = fbh, 
        mc = mc, 
        sbh = sbh, 
        z = z,
        mf_model = 0, 
        sbh_width = 6, 
        nm = nm, 
        Use_S2 = Use_S2, 
        S1_method = S1_method,
        Use_interp = 1,
        S_Tab_Len = S_Tab_Len,
        show_status = 0)
    
    # Print results to a file
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    
    for idx in np.arange(0, nz):
        Str = Str + '  {0:.5E}'.format(r[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()

    PyLab.SaySomething()
    return r

def model_1(theta):
    Use_S2 = 1
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    r = Merger_Rate(
        fbh = fbh, 
        mc = mc, 
        sbh = sbh, 
        z = z,
        mf_model = 0, 
        sbh_width = 6, 
        nm = nm, 
        Use_S2 = Use_S2, 
        S1_method = S1_method,
        Use_interp = 1,
        S_Tab_Len = S_Tab_Len,
        show_status = 0)
    
    # Print results to a file
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    
    for idx in np.arange(0, nz):
        Str = Str + '  {0:.5E}'.format(r[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()

    PyLab.SaySomething()
    return r

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
    
    # Get posterior
    P0 = PyLab.mcmc_derived_stat(model_function = model_0, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    P1 = PyLab.mcmc_derived_stat(model_function = model_1, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    
    np.savez('data/4_merger_rate_posterior.npz', P0 = P0, P1 = P1)
    PyLab.Timer(t1)
