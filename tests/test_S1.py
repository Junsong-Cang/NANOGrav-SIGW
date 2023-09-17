LineWidth = 2
FontSize = 18

reload = 1
fbh = 0.05
mc = 1e-4
mf_model = 0

from src import  merger as mg0
from src import  merger_4 as mg4
from PyLab import *

M = mc + mc

def model_1(sbh):

    S1 = mg4.Get_S1(fbh = fbh, mc = mc, sbh = sbh, M = M, mf_model = mf_model)
    return S1

def model_2(sbh):
    M = mc + mc
    m_ave = mg0.m_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    m2_ave = mg0.m2_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    
    S1 = mg0.S_Factor(fbh = fbh, sbh = sbh, M = M, m_ave = m_ave, m2_ave = m2_ave, t = 3e17, Use_S2=0)
    return S1

n = 100
sbh = np.linspace(0.01, 2, n)
r1 = np.zeros(n)
r2 = np.zeros(n)

if reload:
    
    r1 = Parallel(n_jobs=12)(delayed(model_1)(x) for x in sbh)
    t1 = TimeNow()
    r2 = Parallel(n_jobs=12)(delayed(model_2)(x) for x in sbh)
    Timer(t1)
    np.savez('tmp.npz', r1 = r1, r2 = r2)

r = np.load('tmp.npz')
r1 = r['r1']
r2 = r['r2']
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.plot(sbh, r1, '-k', linewidth=LineWidth, label = 'Analytic')
plt.plot(sbh, r2, '-r', linewidth=LineWidth, label = 'Approximate')

plt.xlabel('$\sigma_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$S_1$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper right')
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

dif = np.sum(1-r2/r1)/n
print(dif)
