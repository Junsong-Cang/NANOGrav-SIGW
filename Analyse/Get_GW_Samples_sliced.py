from src.merger import *

Slice = 3
reload = 0
FileRoot_3 = '/home/dm/gaolq/cjs/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
FileRoot_4 = '/home/dm/gaolq/cjs/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'

v = np.logspace(-3, 10, 80)
ncpu = 54

def Get_params(idx, sample):
    if sample == 3:
        File = FileRoot_3 + '.txt'
    elif sample == 4:
        File = FileRoot_4 + '.txt'
    Tab = np.loadtxt(File)
    chain = Tab[idx,:]
    fbh = 10**chain[7]
    mc = 10**chain[8]
    sbh = chain[9]
    r = fbh, mc, sbh
    return r

def model(idx, sample, Use_S2):
    
    fbh, mc, sbh = Get_params(idx = idx, sample = sample)
    r = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh, 
        v = v,
        mf_model = 0, 
        sbh_width = 7, 
        nm = 50,
        nz = 50,
        show_status = 0,
        Use_interp = 1,
        S1_method = 0,
        Fast = 0,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    # r = np.zeros(len(v))
    if sample == 3:
        PyLab.SaySomething('tmp_3.txt')
    else:
        PyLab.SaySomething('tmp_4.txt')
    return r

if reload:

    n3 = len(np.loadtxt(FileRoot_3 + '.txt'))
    n4 = len(np.loadtxt(FileRoot_4 + '.txt'))
    t1 = PyLab.TimeNow()
    if Slice == 3:
        r30 = Parallel(n_jobs = ncpu)(delayed(model)(idx = idx, sample = 3, Use_S2 = 0) for idx in np.arange(0, n3))
        r31 = Parallel(n_jobs = ncpu)(delayed(model)(idx = idx, sample = 3, Use_S2 = 1) for idx in np.arange(0, n3))
        np.savez('GW_Samples_3.npz', r30 = r30, r31 = r31)
    else:
        r40 = Parallel(n_jobs = ncpu)(delayed(model)(idx = idx, sample = 4, Use_S2 = 0) for idx in np.arange(0, n4))
        r41 = Parallel(n_jobs = ncpu)(delayed(model)(idx = idx, sample = 4, Use_S2 = 1) for idx in np.arange(0, n4))
        np.savez('GW_Samples_4.npz', r40 = r40, r41 = r41)
    PyLab.Timer(t1)
    