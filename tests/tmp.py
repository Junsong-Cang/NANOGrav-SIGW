from p21c_tools import *

nz = 2
z1 = np.linspace(0, 1100, nz)
t = np.zeros(nz)
z2 = np.zeros(nz)

for idx in np.arange(0, nz):
    t[idx] = z2t(z1[idx])

for idx in np.arange(0, nz):
    z2[idx] = t2z(t[idx])
    
print(z1)
print(z2)
