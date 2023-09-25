# compare my merger rate with 2012.02786
from src.merger import *

reload = 1

f = np.logspace(-4, 0, 50)

LineWidth = 2
FontSize = 18

import matplotlib.pyplot as plt
from joblib import Parallel, delayed

f1, r1 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2012.02786.fig1.black.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

def model_1(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 0.6, z = 0, mf_model = 1, sbh_width = 6, nm = 200, Use_S2 = 1, S1_method = 0, Use_interp = 1, S_Tab_Len = 200)
    PyLab.SaySomething()
    return r

f2, r2 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1908.09752.fig2.green_solid.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

def model_2(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 0.02, z = 0, mf_model = 2, sbh_width = 6, nm = 200, Use_S2 = 1, S1_method = 0, Use_interp = 1, S_Tab_Len = 200)
    PyLab.SaySomething()
    return r

f3, r3 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812.01930.fig3.black_solid.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

def model_3(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 0.6, z = 0, mf_model = 0, sbh_width = 6, nm = 200, Use_S2 = 0, S1_method = 0, Use_interp = 1, S_Tab_Len = 200)
    PyLab.SaySomething()
    return r

f4, r4 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812.01930.fig3.black_dashed.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

def model_4(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 2, z = 0, mf_model = 0, sbh_width = 6, nm = 200, Use_S2 = 0, S1_method = 0, Use_interp = 1, S_Tab_Len = 200)
    PyLab.SaySomething()
    return r

f5, r5 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1908.09752.fig2.green_dashed.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

def model_5(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 0.6, z = 0, mf_model = 2, sbh_width = 6, nm = 200, Use_S2 = 0, S1_method = 0, Use_interp = 1, S_Tab_Len = 200)
    PyLab.SaySomething()
    return r

if reload:
    t1 = PyLab.TimeNow()
    r1_ = Parallel(n_jobs = 12)(delayed(model_1)(x) for x in f)
    r2_ = Parallel(n_jobs = 12)(delayed(model_2)(x) for x in f)
    r3_ = Parallel(n_jobs = 12)(delayed(model_3)(x) for x in f)
    r4_ = Parallel(n_jobs = 12)(delayed(model_4)(x) for x in f)
    r5_ = Parallel(n_jobs = 12)(delayed(model_5)(x) for x in f)
    PyLab.Timer(t1)
    np.savez('data/merger_rate.npz', r1_ = r1_, r2_ = r2_, r3_ = r3_, r4_ = r4_, r5_ = r5_)

r = np.load('data/merger_rate.npz')
r1_ = r['r1_']
r2_ = r['r2_']
r3_ = r['r3_']
r4_ = r['r4_']
r5_ = r['r5_']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(f1, r1, '-k', linewidth=LineWidth, label = 'Paper')
plt.loglog(f, r1_, '--k', linewidth=LineWidth, label = 'my codes')

plt.loglog(f2, r2, '-r', linewidth=LineWidth)
plt.loglog(f, r2_, '--r', linewidth=LineWidth)

plt.loglog(f3, r3, '-b', linewidth=LineWidth)
plt.loglog(f, r3_, '--b', linewidth=LineWidth)

plt.loglog(f4, r4, '-g', linewidth=LineWidth)
plt.loglog(f, r4_, '--g', linewidth=LineWidth)

plt.loglog(f5, r5, '-c', linewidth=LineWidth)
plt.loglog(f, r5_, '--c', linewidth=LineWidth)

plt.xlabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[\mathrm{Gpc^{-3} yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
#plt.xlim([-1,7])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
#plt.show()

# print(r2_)
