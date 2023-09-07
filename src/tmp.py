
from merger_module_2 import *

n = 1000
t1 = PyLab.TimeNow()

for idx in np.arange(0, n):
    r = S_factor(fbh = 1e-2, mc = 20, sbh = 0.6, M = 40, nm = 100, nv = 40, sbh_width = 5, method = 0)
PyLab.Timer(t1)

print(r)
r2 = S_factor(fbh = 1e-2, mc = 20, sbh = 0.6, M = 40, nm = 10000, nv = 1000, sbh_width = 10, method = 0)
print(r2)
dif = abs(r - r2)/r2
print(dif)

