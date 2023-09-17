import p21c_tools as p21ct
import scipy.special as sp
from Useful_Numbers import Cosmology as cosmo
from Useful_Numbers import Constants as const
from PyLab import *
try:
    import merger
except:
    from src import merger

Small = 1e-280

def Psi(fbh, mc, sbh, m):
    u = np.log(mc)
    x = np.log(m)
    Chi2 = (x - u)**2 / sbh**2
    Exp = np.exp(-Chi2/2)
    Prefix = fbh/(sbh * m * np.sqrt(2 * np.pi))
    r = Prefix * Exp
    return r

def dR_dm1dm2dm3_kernel(fbh, sbh, mc, m1, m2, m3, m_ave, t):
    
    # Get x
    OmR = 9.1e-5
    zeq = cosmo['zeq']
    zpeq = 1 + zeq
    aeq = 1/zpeq
    msun = const['msun']
    OmM = cosmo['OmM']
    OmL = 1 - OmM - OmR
    G = const['G']
    c = const['c']
    
    RhoEq = cosmo['RhoCr'] * (OmL + OmM*zpeq**3 + OmR * zpeq**4)
    M = m1+m2
    M_SI = M*msun
    x3 = 3 * M_SI /(4 * np.pi * aeq**3 * RhoEq)
    x = x3**(1/3)

    # Get tau
    eta = m1*m2/M**2
    m3_SI = m3*msun

    tau_1 = 384 * aeq**4 * m3**7 * x**4
    tau_2 = 85 * G**3 * eta * M**10
    tau = c**5 * tau_1/tau_2/msun**3
    
    # Get N
    
    N = 0.42 * M/m_ave

    # Get pieces
    G1 = sp.gammainc(58/37, N*(t/tau)**(3/16))
    G2 = sp.gammainc(58/37, N*(t/tau)**(-1/7))
    Gamma = G1 - G2

    p1 = 9 * (t/tau)**(-34/37) / (296 * np.pi * tau)
    p2 = Gamma
    p3 = x**(-3) * N**(53/37) * m_ave**3

    Psi_1 = Psi(fbh = fbh, mc = mc, sbh = sbh, m = m1)
    Psi_2 = Psi(fbh = fbh, mc = mc, sbh = sbh, m = m2)
    Psi_3 = Psi(fbh = fbh, mc = mc, sbh = sbh, m = m3)

    p4 = (Psi_1 * Psi_2 * Psi_3)/(m1 * m2 * m3)

    r = p1 * p2 * p3 * p4
    
    if p3 < -Small:
        raise Exception('kernel is wrong')

    
    return r

def get_m_vec(mc, sbh, sbh_width, nm):
    u = np.log(mc)
    x1 = u - sbh_width*sbh
    x2 = u + sbh_width*sbh
    x = np.linspace(x1, x2, nm)
    m = np.exp(x)
    return m
    
def Get_dR_dx1dx2dx3_3D(fbh, sbh, mc, sbh_width, nm, z, show_status):
    '''
    x is lnm
    '''
    # speed: 8s for nm=100
    t = p21ct.z2t(z)
    m_ave = fbh*mc*np.exp(-sbh**2/2)

    m_vec = get_m_vec(mc = mc, sbh = sbh, sbh_width = sbh_width, nm = nm)
    dR_dx1dx2dx3_vec = np.zeros((nm, nm, nm))
    
    for mid_3 in np.arange(0, nm):
        if show_status:
            print('Get_dR_dx1dx2dx3_vec status:', mid_3/nm)
        for mid_1 in np.arange(0, nm):
            for mid_2 in np.arange(0, nm):

                m1 = m_vec[mid_1]
                m2 = m_vec[mid_2]
                m3 = m_vec[mid_3]

                if mid_1 <= mid_2:
                    
                    dR_dm1dm2dm3 = dR_dm1dm2dm3_kernel(fbh = fbh, sbh = sbh, mc = mc, m1 = m1, m2 = m2, m3 = m3, m_ave = m_ave, t = t)
                    dR_dx1dx2dx3_vec[mid_1, mid_2, mid_3] = m1 * m2 * m3 * dR_dm1dm2dm3
        
        for mid_1 in np.arange(0, nm):
            for mid_2 in np.arange(0, nm):

                if mid_1 > mid_2:
                    dR_dx1dx2dx3_vec[mid_1, mid_2, mid_3] = dR_dx1dx2dx3_vec[mid_2, mid_1, mid_3]

    # r = {'m' : m_vec, 'lnm': np.log(m_vec), 'dR_dx1dx2dx3_3D' : dR_dx1dx2dx3_vec}

    return dR_dx1dx2dx3_vec
                
def Get_dR_dx1dx2dx3_4D(fbh, sbh, mc, sbh_width, nm, zmax, nz, show_status, ncpu = 1):
    
    if nm > 100:
        raise Exception('nm too large, memory might not be enough')

    zp_vec = np.logspace(0, np.log10(zmax + 1), nz)
    z_vec = zp_vec - 1
    
    def Get_3D_Tab(z):
        r = Get_dR_dx1dx2dx3_3D(fbh = fbh, sbh = sbh, mc = mc, sbh_width = sbh_width, nm = nm, z = z, show_status = 0)
        return r
    
    if ncpu == 1:
        dR_dx1dx2dx3_4D = np.zeros((nz, nm, nm, nm))
        
        for zid in np.arange(0, nz):
            z = z_vec[zid]
            if show_status:
                print('Get_dR_dx1dx2dx3_4D status:', zid/nz)
            Tab = Get_3D_Tab(z = z)
            dR_dx1dx2dx3_4D[zid, :,:,:] = Tab[:,:,:]
    else:
        dR_dx1dx2dx3_4D = Parallel(n_jobs=ncpu)(delayed(Get_3D_Tab)(x) for x in z_vec)
    
    m_vec = get_m_vec(mc = mc, sbh = sbh, sbh_width = sbh_width, nm = nm)

    dR_dx1dx2dx3_4D = np.array(dR_dx1dx2dx3_4D)

    r = {'z_vec' : z_vec, 'm_vec' : m_vec, 'dR_dx1dx2dx3_4D' : dR_dx1dx2dx3_4D}

    return r
    # return dR_dx1dx2dx3_4D

def Get_OmGW_kernel(v, Tables, ncpu, show_status):
    
    z_vec = Tables['z_vec']
    m_vec = Tables['m_vec']
    dR_dx1dx2dx3_4D = Tables['dR_dx1dx2dx3_4D']
    
    x_vec = np.log(m_vec)
    nm = len(m_vec)
    nz = len(z_vec)
    
    Is_Scalar(v, 2)
    
    def Find_Fz(zid):
        '''
        Integrate m1, m2, m3
        '''
        
        z = z_vec[zid]
        H = Hubble(z = z)
        zp = 1+z
        vr = zp*v
        dR_dx1dx2dx3 = np.zeros(nm)
        # dR_dx1dx2 = np.zeros(nm)
        fx1x2 = np.zeros(nm)
        fx1 = np.zeros(nm)

        for mid_1 in np.arange(0, nm):
            for mid_2 in np.arange(0, nm):

                for mid_3 in np.arange(0, nm):
                    #print(zid, mid_1, mid_2, mid_3)
                    #print(np.shape(dR_dx1dx2dx3_4D))
                    dR_dx1dx2dx3[mid_3] = dR_dx1dx2dx3_4D[zid, mid_1, mid_2, mid_3]
                    #print(dR_dx1dx2dx3[mid_3])
                
                dR_dx1dx2 = np.trapz(x = x_vec, y = dR_dx1dx2dx3)
                if True in (dR_dx1dx2dx3 < -Small):
                    raise Exception('dR_dx1dx2dx3 is wrong')
                m1 = m_vec[mid_1]
                m2 = m_vec[mid_2]
                M = m1 + m2
                eta = m1*m2/M**2
                dEGW_dvr = merger.Get_dEGW_dv(v = vr, eta = eta, M = M)
                
                fx1x2[mid_2] = dR_dx1dx2 * dEGW_dvr / (zp*H)
                if dR_dx1dx2 < -Small:
                    raise Exception('dR_dx1dx2 is wrong')
        
                if fx1x2[mid_2] < -Small:
                    raise Exception('fx1fx2 is wrong')
        
            
            fx1[mid_1] = np.trapz(x = x_vec, y = fx1x2)
            if True in (fx1 < -Small):
                raise Exception('fx1 is wrong')
        
        fx = np.trapz(x = x_vec, y = fx1)
        
        # Need some prefactors
        RhoCrC2 = 7.730937688449169e-10
        Prefix = v/RhoCrC2
        
        r = Prefix * fx
        
        if fx < -Small:
            raise Exception('r is wrong')
        
        return r
    
    if ncpu == 1:
        fz = np.zeros(nz)
        for zid in np.arange(0, nz):
            if show_status:
                print('Get_OmGW_kernel status:', zid/nz)
            fz[zid] = Find_Fz(zid)
    else:
        idxs = np.arange(0, nz)
        fz = Parallel(n_jobs = ncpu)(delayed(Find_Fz)(idx) for idx in idxs)
    
    fz = np.array(fz)
    zp = z_vec + 1
    OmGW = np.trapz(x = np.log(zp), y = zp*fz)
    
    if True in (fz < -Small):
        raise Exception('fz is wrong')
    
    if (np.min(zp) < -Small):
        raise Exception('zp is wrong')
    
    if OmGW < -Small:
        raise Exception('OmGW is negative')
    
    return OmGW
    
def Get_OmGW(
        fbh = 0.002930716205753681,
        mc = 30,
        sbh = 1,
        v = np.logspace(-6, 4, 100),
        sbh_width = 5,
        zmax = 50,
        nm = 50,
        nz = 100,
        show_status = 0,
        ncpu = 12):
    
    if show_status:
            print('Get_OmGW : calling Get_dR_dx1dx2dx3_4D')
    Tables = Get_dR_dx1dx2dx3_4D(fbh = fbh, sbh = sbh, mc = mc, sbh_width = sbh_width, 
                                 nm = nm, zmax = zmax, nz = nz, show_status = show_status, ncpu = ncpu)
    # return Tables

    if show_status:
            print('Get_OmGW : Get_dR_dx1dx2dx3_4D call complete')
    
    nv = len(v)
    r = np.zeros(nv)
    
    for idx in np.arange(0, nv):

        if show_status:
            print('Get_OmGW status:', idx/nv)
        
        v_ = v[idx]
        r[idx] = Get_OmGW_kernel(v = v_, Tables = Tables, ncpu = ncpu, show_status = 0)
    
    return r
