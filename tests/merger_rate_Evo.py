# compare my merger rate evolution with 2012.02786 and 1812

from src.merger import *

LineWidth = 2
FontSize = 18
nz = 1000

import matplotlib.pyplot as plt


z1, r1 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812_01930_fig10_red.txt',
    nx = nz,
    model = 1,
    Convert_x = 0,
    Convert_y = 1)

zp2, r2 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2012_02786_fig6a_black.txt',
    nx = nz,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

z2 = zp2-1
zp1 = 1+z1

r3 = np.zeros(nz)
for idx in np.arange(0, nz):
    r3[idx] = Merger_Rate_Scale_Factor(z = z2[idx], fbh = 1e-2, Use_S2 = 0)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(zp1, r1, '-k', linewidth=LineWidth, label = '1812.01930')
plt.loglog(zp2, r2, '-r', linewidth=LineWidth, label = '2012.02786')
plt.loglog(zp2, r3, '-b', linewidth=LineWidth, label = 'My codes')

plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$R(z)/R(0)$',fontsize=FontSize,fontname='Times New Roman')
plt.title('Merger Rate Scale Factor',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
#plt.xlim([-1,7])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
