LineWidth = 2
FontSize = 18

from PyLab import *

f1, r1 = Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812.01930.fig3.black_solid.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f2, r2 = Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812.01930.fig3.black_dashed.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(f1, r1, '-k', linewidth=LineWidth, label = '$\sigma = 0.6$')
plt.loglog(f2, r2, '--k', linewidth=LineWidth, label = '$\sigma = 2$')

plt.xlabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[\mathrm{Gpc^{-3} yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

a = np.zeros((100,2))
a[:,0] = f2[:]
a[:,1] = r2[:]
file = '/Users/cangtao/Desktop/tmp.txt'

np.savetxt(fname = file, X = a, fmt = '%0.4e', delimiter = '    ')
