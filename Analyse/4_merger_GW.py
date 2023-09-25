from src.merger import *

reload = 1
v = np.logspace(-3, 10, 100)

LineWidth = 2
FontSize = 18

params = np.array(
    [[0.1, 10**-3.589, 0.54],
    [0.1, 10**-3.589, 0.01],
    [0.1, 10**-3.589, 1],
    [0.1, 10**-3.589, 2]])

import matplotlib.pyplot as plt
n = len(params)

def model(f, m, s):
    # r = Get_dOmGW_dlnv(fbh = f, mc = m, sbh = s, v = v, mf_model = 0, sbh_width = 6, nm = 100, nz = 50, ncpu = 12)
    if s < 0.05:
         mf_model = 2
    else:
         mf_model = 0
    r = Get_dOmGW_dlnv(
        fbh = f, 
        mc = m,
        sbh = s, 
        v = v,
        mf_model = mf_model,
        sbh_width = 7,
        nm = 50,
        nz = 100,
        show_status = 0,
        Use_interp = 1,
        S1_method = 0,
        Fast = 0,
        S_Tab_Len = 200,
        Use_S2 = 0,
        Precision = 0,
        ncpu = 11)
    return r

if reload:
    
    r = np.zeros((n, len(v)))

    t1 = PyLab.TimeNow()
    for idx in np.arange(0, n):
        print('status : ', idx/n)
        fbh = params[idx, 0]
        mc = params[idx, 1]
        sbh = params[idx, 2]
        r[idx, :] = model(f = fbh, m = mc, s = sbh)
    np.savez('data/4_merger_GW.npz', r = r, v = v)
    PyLab.Timer(t1)

r = np.load('data/4_merger_GW.npz')['r']

print('Parameter Values:')
print('    fbh        mc        sigma')
for idx in np.arange(0, n):
    fbh = params[idx, 0]
    mc = params[idx, 1]
    sbh = params[idx, 2]
    print("{0:.4E}".format(fbh), " {0:.4E}".format(mc), "  {0:.4f}".format(sbh))
