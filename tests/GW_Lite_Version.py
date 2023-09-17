from src.merger import *
from PyLab import *

m = 1e-2
fbh = 1e-1
v = np.logspace(0, 8, 40)
Use_S2 = 1
S1_method = 0

t1 = TimeNow()
g1 = Get_dOmGW_dlnv(
    fbh = fbh, 
    mc = m, 
    v = v, 
    mf_model = 2, 
    nz = 10000,
    zmax = 100000, 
    S1_method = S1_method, 
    Use_S2 = Use_S2)
Timer(t1)
t1 = TimeNow()
g2 = dOmGW_dlnv_Lite_kernel_pro(
    fbh = fbh,
    mc = m,
    v_vec =v,
    Use_S2 = Use_S2,
    S1_method = S1_method,
    nz = 50,
    Precision = 1e-3,
    ncpu = 12,
    show_status=0)
Timer(t1)

dif = np.sum(np.abs(g1 - g2))/np.sum(g1)/len(v)
print(dif)

plt.loglog(v, g1, 'k')
plt.loglog(v, g2, '--r')
plt.show()

