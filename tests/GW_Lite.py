from src.merger import *
from PyLab import *

'''
n = 100000
# v = np.logspace(-9, np.log10(11572.365475258206), n)
v = np.logspace(-20, 6, n)
z = 100
m = 1

vr = (1+z)*v
dEdvr = np.zeros(n)
dEdv = np.zeros(n)

for idx in np.arange(0,n):
    v_ = v[idx]
    vr_ = v[idx] * (1+z)
    dEdvr[idx] = Get_dEGW_dv(v = vr_, eta = 1/4, M = 2*m)
    
    dEdv[idx] = Get_dEGW_dv(v = v_, eta = 1/4, M = 2*m)


E = \int dvr dE/dvr
E = \int dv dE/dv

E1 = np.trapz(x = np.log(vr), y = vr * dEdvr)
E2 = np.trapz(x = np.log(v), y = v * dEdv)

E1 = np.trapz(x = vr, y = dEdvr)
E2 = np.trapz(x = v, y = dEdv)

print(E1)
print(E2)
'''

n = 100
v = np.linspace(8100.451859543111, 11572.365475258206, n)
fbh = 0.05
m = 1

vv = 2856
z = 2
rr = Get_dEGW_dv(v = vv * (1+z), eta = 1/4, M = 2*m)
print(rr)

'''
r = np.zeros(n)
for idx in np.arange(0, n):
    r[idx] = Get_dEGW_dv(v = v[idx], eta = 1/4, M = 2*m)
plt.loglog(v, r)
# plt.show()
'''

vvv = np.array([1e2, 8000])
rrr = Get_dOmGW_dlnv(fbh = fbh, mc = m, v = vvv, mf_model=2, nz = 10000, zmax = 1000)
print(rrr)
