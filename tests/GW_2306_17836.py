
Use_S2 = 0
S1_method = 0
LineWidth = 2
FontSize = 16

reload = 1
ResultFile='/Users/cangtao/Desktop/tmp.pdf'

from src.merger import *
from PyLab import *

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
    r1 = Get_dOmGW_dlnv(
        fbh = 1,
        mc = 1e5,
        v = v,
        mf_model = 2,
        nz = 50,
        show_status = 1,
        S1_method = S1_method,
        Use_S2 = Use_S2,
        Precision = 0)

    np.savez('tmp.npz', r1 = r1)

r1 = np.load('tmp.npz')['r1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(v, r0, '-k', linewidth=LineWidth, label = '2306.17836')
plt.loglog(v, r1, '-r', linewidth=LineWidth, label = 'My Codes')

plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('${\mathrm{d}}\Omega_{\mathrm{GW}}/{\mathrm{dln}}f$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper left')
plt.tight_layout()

plt.savefig(ResultFile)
print('Plot saved to :')
print(ResultFile)
