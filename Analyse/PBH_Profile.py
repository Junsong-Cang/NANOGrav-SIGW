from src.main import *

LineWidth = 2
FontSize = 18
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Phi_compare.pdf'

p1 = Phi_Profile(
    A = 10**-0.81,
    Sigma = 0.02,
    kc = 10**8.1)

p2 = Phi_Profile(
    A = 10**-0.81,
    Sigma = 0.5,
    kc = 10**8.1)

p3 = Phi_Profile(
    A = 10**-0.81,
    Sigma = 1.25,
    kc = 10**8.1)

def Normalise_Profile(P):
    x = P['m_vec']/P['mc']
    ya = P['Phi_vec']/P['fbh']
    yb = P['Best_Fit']/P['fbh']
    return x, ya, yb

x1, ya1, yb1 = Normalise_Profile(p1)
x2, ya2, yb2 = Normalise_Profile(p2)
x3, ya3, yb3 = Normalise_Profile(p3)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.semilogx(x3, ya3, '-k', linewidth=LineWidth, label = '$\Delta = 1.25$')
plt.semilogx(x2, ya2, '-r', linewidth=LineWidth, label = '$\Delta = 0.5$')
plt.semilogx(x1, ya1, '-b', linewidth=LineWidth, label = '$\Delta = 0.02$')

plt.semilogx(x3, yb3, '--k', linewidth=LineWidth)
plt.semilogx(x2, yb2, '--r', linewidth=LineWidth)
plt.semilogx(x1, yb1, '--b', linewidth=LineWidth)

plt.xlabel('$m/m_{\mathrm{c}}$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\psi/f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

#plt.title('A plot example',fontsize=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper left')
#plt.xlim([1e-1,7])
plt.ylim([0,1.6])
plt.tight_layout()

plt.savefig(ResultFile, dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')
print('Plot saved to :')
print(ResultFile)
