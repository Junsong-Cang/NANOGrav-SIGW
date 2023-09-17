from src.merger import *

reload = 1
v = np.logspace(0, 10, 100)

LineWidth = 2
FontSize = 18

fbh = np.array([1e-2, 1e-3, 0.1])
mc = np.array([1e-3, 1e-4, 1e-2])
sbh = np.array([0.3, 0.1, 1])

import matplotlib.pyplot as plt
def model(f, m, s):
    # r = Get_dOmGW_dlnv(fbh = f, mc = m, sbh = s, v = v, mf_model = 0, sbh_width = 6, nm = 100, nz = 50, ncpu = 12)
    r = Get_dOmGW_dlnv(
        fbh = f, 
        mc = m,
        sbh = s, 
        v = v,
        mf_model = 0,
        sbh_width = 7,
        nm = 50,
        nz = 50,
        show_status = 1,
        Use_interp = 1,
        S1_method = 0,
        Fast = 0,
        S_Tab_Len = 200,
        Use_S2 = 1,
        Precision = 1e-2,
        ncpu = 12)
    return r

if reload:
    
    r = np.zeros((9, len(v)))

    idx = 0
    t1 = PyLab.TimeNow()
    r[0, :] = model(f = fbh[0], m = mc[0], s = sbh[0])

    for fid in np.arange(1,3):
        for mid in np.arange(1,3):
            for sid in np.arange(1,3):
                print('--------', idx,'--------')
                PyLab.SaySomething()
                idx = idx + 1
                f_ = fbh[fid]
                m_ = mc[mid]
                s_ = sbh[sid]
                r[idx, :] = model(f = f_, m = m_, s = s_)
    np.savez('data/3_merger_lite.npz', r = r)
    PyLab.Timer(t1)

r = np.load('data/3_merger_lite.npz')['r']

print('Parameter Values:')
print('    fbh        mc        sigma')

print("{0:.4E}".format(fbh[0]), "{0:.4E}".format(mc[0]), "  {0:.4f}".format(sbh[0]))

for fid in np.arange(1,3):
        for mid in np.arange(1,3):
            for sid in np.arange(1,3):
                f = fbh[fid]
                m = mc[mid]
                s = sbh[sid]
                print("{0:.4E}".format(f), "{0:.4E}".format(m), "  {0:.4f}".format(s))
                
    
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()
fig.set_size_inches(8, 4)

# plt.loglog(v, r[0,:], 'k', linewidth=LineWidth, label = '$f_{\mathrm{bh}} = 10^{-2},\ m_{\mathrm{c}} = 10^{-3}m_{\odot},\ \sigma_{\mathrm{bh}} = 0.3$')
plt.loglog(v, r[0,:], 'k', linewidth=LineWidth, label = '$[10^{-2},\ 10^{-3},\ 0.3]$')
plt.loglog(v, r[1,:], 'r', linewidth=LineWidth, label = '$[10^{-3},\ 10^{-4},\ 0.1]$')
plt.loglog(v, r[2,:], 'b', linewidth=LineWidth, label = '$[10^{-3},\ 10^{-4},\ 1]$')
plt.loglog(v, r[3,:], 'c', linewidth=LineWidth, label = '$[10^{-3},\ 10^{-2},\ 0.1]$')
plt.loglog(v, r[4,:], 'g', linewidth=LineWidth, label = '$[10^{-3},\ 10^{-3},\ 1]$')
plt.loglog(v, r[5,:], 'y', linewidth=LineWidth, label = '$[0.1,\ 10^{-4},\ 0.1]$')
plt.loglog(v, r[6,:], 'm', linewidth=LineWidth, label = '$[0.1,\ 10^{-4},\ 1]$')
plt.loglog(v, r[7,:], 'brown', linewidth=LineWidth, label = '$[0.1,\ 10^{-2},\ 0.1]$')
plt.loglog(v, r[8,:], 'purple', linewidth=LineWidth, label = '$[0.1,\ 10^{-2},\ 1]$')

plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize = FontSize,loc='center left', bbox_to_anchor=(1, 0.5))

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

#plt.xlim([-1,7])
#plt.ylim([-1.2,1.2])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

gg = np.interp(x = 1, xp = np.log10(v), fp = r[0,:])
print(gg)
print(np.shape(r))
n = len(v)
file = '/Users/cangtao/Desktop/tmp.txt'

d = np.zeros((10, n))
d[0,:] = v[:]
d[1:10, :] = r[:,:]

np.savetxt(fname = file, X = d, fmt = '%.8E', delimiter = '    ')
