
Use_S2 = 0
S1_method = 0
LineWidth = 2
FontSize = 18

reload = 1

from src import merger as mg
from PyLab import *
from src.merger_2306 import *

# Load paper results
F = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2306_17836_fig1b_dotted.txt'
v, r0 = Read_Curve(
    File = F,
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

h2 = 0.6766**2
r0 = r0/h2

if reload:
    r1 = mg.Get_dOmGW_dlnv(
        fbh = 1,
        mc = 1e5,
        v = v,
        mf_model = 2,
        nz = 200,
        show_status = 1,
        S1_method = S1_method,
        Use_S2 = Use_S2,
        Precision = 0)

    np.savez('tmp.npz', r1 = r1)
    
r= np.load('tmp.npz')
r1 = r['r1']

t1 = TimeNow()
r2 = Get_OmGW(v)
Timer(t1)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(v, r0, '-k', linewidth=LineWidth, label = 'Paper')
plt.loglog(v, r1, '-r', linewidth=LineWidth, label = 'My codes')
#plt.loglog(v, r2, '-b', linewidth=LineWidth, label = 'Re-produced')
plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
#plt.xlim([-1,7])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

#print(r0)
