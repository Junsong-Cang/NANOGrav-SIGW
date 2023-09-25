# compare my merger rate with 2012.02786
from src.merger import *

reload = 1

LineWidth = 2
FontSize = 18

Use_S2 = 0
S1_method = 0
fbh = 5e-3
fid = 5
m = 1e12
scale = 1

nz = 100

import matplotlib.pyplot as plt
h2 = 0.6766**2

files = [
    '2307_01457_fig2_green.txt', # 0
    '2307_01457_fig2_dark_blue.txt', # 1
    '2307_01457_fig2_bright_blue.txt', # 2
    '2307_01457_fig2_yellow.txt', # 3
    '2307_01457_fig2_brown.txt', # 4
    '2307_01457_fig2_red.txt'] # 5

v1, g1 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/' + files[fid],
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)
g1_ = Get_dOmGW_dlnv(
        fbh = fbh, 
        mc = m,
        sbh = 1.0, 
        v = v1,
        mf_model = 2, 
        nz = 100,
        show_status = 0,
        Use_interp = 1,
        S1_method = 0,
        Fast = 0,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 12)

print(g1_)
g1_ = g1_*h2

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(v1, g1, 'k', linewidth=LineWidth, label = '2307.01457')
plt.loglog(v1, g1_/scale, '--k', linewidth=LineWidth, label = 'my code')

plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}h^2$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper left')
#plt.xlim([1, 1e6])
#plt.ylim([1e-12, 1e-6])

plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
