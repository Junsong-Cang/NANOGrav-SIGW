from src.merger_1707 import *

reload = 0
fbh = 0.002930716205753681
mc = 30
sbh = 1
v = np.logspace(-6, 4, 100)
ncpu = 12

if reload:
    t1 = TimeNow()
    g = Get_OmGW(fbh = fbh, mc = mc, sbh = sbh, sbh_width = 6, zmax = 50, nm = 30, nz = 100, show_status = 1, ncpu = ncpu)
    Timer(t1)
    np.savez('tmp.npz', g = g)

r = np.load('tmp.npz')

v2, g2 = Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1707.01480.fig2.red_solid_lower.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

g = r['g']
g = -g * 8e7
#g = -g * 31557600
plt.loglog(v2, g2, 'k')
plt.loglog(v, g, 'r')

print(g)
plt.show()
