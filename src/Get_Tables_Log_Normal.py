
reload = 1
DataFile = 'data/GW_Tab_logNormal.npz'
ncpu = 12
nf = 20
nfc = 30
ns = 40

nu = 500
nv = 100

from src.src_2 import *
f = np.logspace(-8.6, -7.6, nf)
fc = np.logspace(-9.6, -6.6, nfc)
s = np.logspace(-2, 1, ns)

from joblib import Parallel, delayed

c = 1.546e-15
k = f/c
kc = fc/c

params = []
for kid in np.arange(0, nfc):
    for sid in np.arange(0, ns):
        p = [kc[kid], s[sid]]
        params.append(p)

def Get_Sample(idx = 2):
    theta = params[idx]
    kc_, Sigma = theta
    os.system('echo -- >> tmp.txt')
    r = np.linspace(0, 1, nf)
    for kid in np.arange(0, nf):
        r[kid] = dOmGW_dlnk(A=1.0, kc = kc_, Sigma = Sigma, k = k[kid], nu = nu, nv = nv)
    return r

if reload:
    t1 = time.time()
    xs = np.arange(0, ns*nfc)
    fk = Parallel(n_jobs = ncpu)(delayed(Get_Sample)(x) for x in xs)
    t2 = time.time()
    print('Time = ', t2 - t1)
    np.savez(DataFile, fk = fk)
else:
    pass

fk = np.load(DataFile)['fk']
s = np.shape(fk)

print(s)
