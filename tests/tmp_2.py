from src.merger import *
from PyLab import *

n = 1000

z = np.logspace(0, 5, n) - 1
r = np.zeros(n)
t = np.zeros(n)

for idx in np.arange(0,n):
    r[idx] = Merger_Rate_Lite(fbh = 1, mc = 1e5, z = z[idx], Use_S2 = 0, S1_method = 1)
    t[idx] = p21ct.z2t(z = z[idx], unit=1)

plt.loglog(t, r)
plt.show()
