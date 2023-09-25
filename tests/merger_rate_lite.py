from src.merger import *

reload = 1
sbh = 0.5
fbh = 1e-5
nm = 100
mc = np.logspace(-7, 4, nm)

def model(m, mf):
    r = Merger_Rate(
        fbh = fbh, 
        mc = m,
        sbh = sbh,
        z = 0, 
        mf_model = mf, 
        sbh_width = 6, 
        nm = 50,
        Use_S2 = 0, 
        S1_method = 0,
        Use_interp = 1,
        S_Tab_Len = 200,
        show_status = 0)
    return r

if reload:
    t1 = PyLab.TimeNow()
    r1 = Parallel(n_jobs = 12)(delayed(model)(m = x, mf = 0) for x in mc)
    r2 = Parallel(n_jobs = 12)(delayed(model)(m = x, mf = 2) for x in mc)
    PyLab.Timer(t1)
    np.savez('tmp.npz', r1 = r1, r2 = r2)
