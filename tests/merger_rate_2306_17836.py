
Use_S2 = 0
S1_method = 0
LineWidth = 2
FontSize = 18

reload = 1

from src import merger as mg
from PyLab import *

if reload:
    
    # Load paper results
    F = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2306_17836_fig1a_dotted.txt'
    t, r0 = Read_Curve(
        File = F,
        nx = 200,
        model = 2,
        Convert_x = 1,
        Convert_y = 1)
    yr = 365.25 * 24 * 3600
    t = t * yr
    n = len(t)
    
    z = np.zeros(n)
    r1 = np.zeros(n)
    
    # Get my results
    for idx in np.arange(0, n):
        print(idx/n)
        t_ = t[idx]
        z_ = mg.p21ct.t2z(t = t[idx])
        z[idx] = z_
        r1[idx] = mg.Merger_Rate_Lite(fbh = 1, mc = 1e5, z = z[idx], Use_S2 = Use_S2, S1_method = S1_method)
    
    np.savez('tmp.npz', z = z, r0 = r0, r1 = r1)
    
r= np.load('tmp.npz')
z = r['z']
r0 = r['r0']
r1 = r['r1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(1+z, r0, '-k', linewidth=LineWidth, label = 'Paper')
plt.loglog(1+z, r1, '-r', linewidth=LineWidth, label = 'My codes')

plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[\mathrm{Gpc^{-3} yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
#plt.xlim([-1,7])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

#print(r0)
print(z.min())
