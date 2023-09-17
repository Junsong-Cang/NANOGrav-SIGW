from src.merger import *

fbh = 1e-2
mc = 1e-3
sbh = 0.3
v = np.logspace(0, 10, 100)
LineWidth = 2
FontSize = 18

reload = 0

from PyLab import *

if reload:
    # Time is about 250s
    t1 = TimeNow()
    r1 = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh,
        v = v,
        mf_model = 0, 
        sbh_width = 6, 
        nm = 100,
        nz = 1000,
        zmax = 10**12,
        show_status = 1,
        Use_interp = 1,
        S1_method = 0,
        Fast = 1,
        S_Tab_Len = 200,
        Use_S2 = 1,
        Precision = 1e-2,
        ncpu = 12)
    Timer(t1)
    np.savez('tmp.npz', r1 = r1)

r1 = np.load('tmp.npz')['r1']

t1 = TimeNow()
r2 = Get_dOmGW_dlnv(
    fbh = fbh,
    mc = mc,
    sbh = sbh,
    v = v,
    mf_model = 0, 
    sbh_width = 6, 
    nm = 50,
    nz = 50,
    zmax = 50,
    show_status = 1,
    Use_interp = 1,
    S1_method = 0,
    Fast = 1,
    S_Tab_Len = 200,
    Use_S2 = 1,
    Precision = 1e-3,
    ncpu = 12)
Timer(t1)

t1 = TimeNow()
r3 = Get_dOmGW_dlnv(
    fbh = fbh,
    mc = mc,
    sbh = sbh, 
    v = v,
    mf_model = 0, 
    sbh_width = 6, 
    nm = 100,
    nz = 100,
    zmax = 1000,
    show_status = 1,
    Use_interp = 1,
    S1_method = 0,
    Fast = 0,
    S_Tab_Len = 200,
    Use_S2 = 1,
    Precision = 0,
    ncpu = 12)
Timer(t1)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(v, r1, '-k', linewidth=LineWidth, label = 'HiRes')
plt.loglog(v, r2, '-r', linewidth=LineWidth, label = 'Fast')
plt.loglog(v, r3, '+b', linewidth=LineWidth, label = 'Pro')

plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'lower left')
#plt.xlim([-1,7])
#plt.ylim([-1.2,1.2])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
