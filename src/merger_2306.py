from PyLab import *
try:
    import merger as mg
except:
    import src.merger as mg

# Load paper results
def Load_Merger_Rate():
    F = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2306_17836_fig1a_dotted.txt'
    t, R = Read_Curve(
        File = F,
        nx = 10000,
        model = 2,
        Convert_x = 1,
        Convert_y = 1)
    yr = 365.25 * 24 * 3600
    t = t * yr
    nz = len(t)    
    z = np.zeros(nz)

    # Get my results
    for idx in np.arange(0, nz):
        z[idx] = mg.p21ct.t2z(t = t[idx])
    return z, R

z_vec, R_vec = Load_Merger_Rate()

def Merger_Rate(z):
    x = np.log(1+z)
    xp = np.log(1+z_vec)
    fp = R_vec
    R = np.interp(x = x, xp = xp, fp = fp)
    return R

def Get_OmGW_kernel(v):
    
    Gpc3yr = 9.274526196872258e+83
    nz = len(z_vec)
    vr_vec = v*z_vec
    zp_vec = z_vec + 1
    H_vec = Hubble(z_vec)
    dEGW_dvr = np.zeros(nz)

    for idx in np.arange(0, nz):
        zp = z_vec[idx] + 1
        vr = vr_vec[idx]
        dEGW_dv = mg.Get_dEGW_dv(v = vr, eta = 1/4, M = 2e5)
        dEGW_dvr[idx] = dEGW_dv/zp
    
    R = R_vec / Gpc3yr

    RhoCrC2 = 7.730937688449169e-10
    p1 = v/RhoCrC2
    p2 = R * dEGW_dvr/(H_vec * (1+z_vec))
    dOmGW_dzp = p1*p2
    
    # Ok z is in decrasing order
    x = np.log(zp_vec[::-1])
    y = zp_vec * dOmGW_dzp
    y = y[::-1]

    r = np.trapz(x = x, y = y)
    
    return r

def Get_OmGW(v):
    
    nv = len(v)
    r = np.zeros(nv)
    
    for idx in np.arange(0, nv):
        r[idx] = Get_OmGW_kernel(v = v[idx])
    
    return r

