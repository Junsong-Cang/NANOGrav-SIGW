import numpy as np
import scipy.special as sp
import p21c_tools as p21ct
import PyLab
import warnings
import scipy.io
from joblib import Parallel, delayed
from Useful_Numbers import Cosmology as cosmo
from Useful_Numbers import Constants as const
from src.merger_lite_debug import *

# load kummerU template computed by matlab, I donno why but python result is strange!
# a simple script like this can generate the template:
'''
clear
fbh_vec = logspace(-20, 0, 10000);
sigma_M = 0.004;
for idx=1:length(fbh_vec)
    U_vec(idx) = kummerU(21/74, 1/2, 5*fbh_vec(idx)^2/(6*sigma_M^2));
end
save /Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template.mat fbh_vec U_vec sigma_M
'''

kummerU_template = scipy.io.loadmat('/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template.mat')

print('--------dEdvr has (1+z) term now, check again!!!---------')

def Psi(mc, sbh, m, mf_model):
    '''
    dnbh/dlnm/nbh
    ----inputs----
    mf_model : MF model
        0 : log-normal in rho
        1 : log-normal in nbh
        2 : monochromatic
    '''
    if mf_model == 1:
        x = np.log(m)
        u = np.log(mc)
        chi2 = ((x-u)/sbh)**2
        result = np.exp(-chi2/2)/(sbh*np.sqrt(2 * np.pi))
    elif mf_model == 0:
        x = np.log(m)
        u = np.log(mc)
        chi2 = ((x-u)/sbh)**2
        Phi = np.exp(-chi2/2)/(sbh*np.sqrt(2 * np.pi))
        # C is the result given by mathematica, best do some numerical checks
        C = mc/(m*np.exp(sbh**2/2))
        result = Phi * C
    return result

def OmGW_2_h(OmGW, f):
    H0 = 2.192695336552484e-18
    C = 4 * np.pi**2 /(3 * H0**2)
    h2 = OmGW/C/f**2
    if True in h2<0:
        raise Exception
    h = h2**0.5
    return h

def h2_OmGW(h, f):
    H0 = 2.192695336552484e-18
    C = 4 * np.pi**2 /(3 * H0**2)
    OmGW = C * f**2 * h**2
    return OmGW

def m_ave_fun(mc, sbh, mf_model):
    '''
    Mean of m
    1 cpu speed: 
        20000 calls / s (mf_model = 1)
        2000000 calls / s (mf_model = 0)
    '''
    sbh_width = 10
    nm = 2000
    
    if mf_model == 0:
        u = np.log(mc)
        x1 = u - sbh_width*sbh
        x2 = u + sbh_width*sbh
        x = np.linspace(x1, x2, nm)
        m = np.exp(x)
        mf = Psi(mc = mc, sbh = sbh, m = m, mf_model = 0)
        norm = np.trapz(x = x, y = mf)
        if abs(norm-1) > 0.1:
            raise Exception('Mass function is not convergent, please increase precision by setting sbh_width and nm')
        fun = m*mf
        result = np.trapz(x = x, y = fun)
    elif mf_model == 1:
        result = mc*np.exp(sbh**2/2)
    return result

def m2_ave_fun(mc, sbh, mf_model):
    '''
    Mean of m^2
    '''
    sbh_width = 10
    nm = 2000

    if mf_model == 0:
        u = np.log(mc)
        x1 = u - sbh_width*sbh
        x2 = u + sbh_width*sbh
        x = np.linspace(x1, x2, nm)
        m = np.exp(x)
        mf = Psi(mc = mc, sbh = sbh, m = m, mf_model = 0)
        norm = np.trapz(x = x, y = mf)
        if abs(norm-1) > 0.1:
            raise Exception('Mass function is not convergent, please increase precision by setting sbh_width and nm')
        fun = m**2 * mf
        result = np.trapz(x = x, y = fun)
    elif mf_model == 1:
        result = mc**2 * np.exp(2 * sbh**2)
    return result

def Get_S1_analytic(fbh, m_ave, mc, sbh, M, nm, nv, mf_model):
    
    # Speed: 30 calls / s
    sbh_width_s1 = 6
    # S1_nm = 200
    # S1_nv = 100

    sigma_M = kummerU_template['sigma_M'][0][0]
    
    Ny = M*fbh/(m_ave * (fbh + sigma_M))
    
    def int_m_kernel(v):
        if mf_model == 2:
            z4f = mc*v/(m_ave*Ny)
            F = sp.hyp2f1(-0.5, 0.75, 1.25, -9*z4f**2/16) - 1
            r = F
            return r
        u = np.log(mc)
        x1 = u - sbh_width_s1 * sbh
        x2 = u + sbh_width_s1 * sbh
        x_vec = np.linspace(x1, x2, nm)
        m_vec = np.exp(x_vec)
        mf = Psi(mc = mc, sbh = sbh, m = m_vec, mf_model = mf_model)
        z4f = m_vec*v/(m_ave*Ny)
        F = sp.hyp2f1(-0.5, 0.75, 1.25, -9*z4f**2/16) - 1
        fun = mf * F
        r = np.trapz(x = x_vec, y = fun)
        return r
    
    def int_m_vec(v):
        n = len(v)
        r = np.zeros(n)
        for idx in np.arange(0, n):
            r[idx] = int_m_kernel(v[idx])
        return r
    
    def int_v_kernel(v):
        p1 = np.exp(-Ny)/sp.gamma(21/37)
        p2 = v**(-16/37)
        
        p3 = Ny*int_m_vec(v)
        p4 = 0.3 * (sigma_M * v/fbh)**2

        p5 = np.exp(-p3-p4)
        
        r = p1 * p2 *p5
        
        return r
    
    v, fv = PyLab.Map(
        F = int_v_kernel,
        Start = 0,
        Width = 1,
        MinX = -np.inf,
        MaxX = np.inf,
        nx = nv,
        Precision = 5e-3,
        Max_Iteration = 50,
        Use_log_x = 0,
        flat_count = 20)
    
    r = np.trapz(x = np.log(v), y = v * fv)
    
    return r

def Get_S1_approximate(fbh, sbh, M, m_ave, m2_ave):
    
    # Speed: 3 * 10^4 calls / s

    sigma_M = kummerU_template['sigma_M'][0][0]
    Ny = M*fbh/(m_ave*(fbh + sigma_M))
    
    # Get C
    if sbh >2:
        # raise Exception('Approximation for C will be compromised for sbh > 2!')
        warnings.warn('Approximation for C will be compromised for sbh > 2!')
    c1 = (fbh/(sigma_M*m_ave))**2 * m2_ave
    c2a = sp.gamma(29/37)/np.sqrt(np.pi)
    log_fbh_axis = np.log(kummerU_template['fbh_vec'][0,:])
    U_axis = kummerU_template['U_vec'][0,:]
    U = np.interp(x = np.log(fbh), xp = log_fbh_axis, fp = U_axis)
    # U = sp.hyp1f1(21/74, 1/2, 5*fbh**2/(6*sigma_M**2))
    
    c2 = (c2a*U)**(-74/21)-1
    C = c1/c2

    # Get S1
    s1a = 1.42*np.exp(-Ny)
    s1b = m2_ave/(m_ave**2)/(Ny + C)
    s1c = (sigma_M/fbh)**2
    S1 = s1a * (s1b + s1c)**(-21/74)

    return S1

def S_Factor(fbh, sbh, mc, M, m_ave, m2_ave, mf_model, t, Use_S2, S1_method):
    '''
    Sacrificed simplicity for the sake of maximum speed
    ----inputs----
    S1_method : 0 - analytic
                1 - approximate
                2 - no S1
    '''
    
    if S1_method == 0:
        S1 = Get_S1_analytic(fbh = fbh, mc = mc, sbh = sbh, M = M, mf_model = mf_model, m_ave = m_ave, nm = 200, nv = 100)
    elif S1_method == 1:
        S1 = Get_S1_approximate(fbh = fbh, sbh = sbh, M = M, m_ave = m_ave, m2_ave = m2_ave)
    else:
        S1 = 1
    
    if Use_S2:
        t0 = 4.3523939036782445e+17
        fbh_ = fbh*(t/t0)**0.44
        S2= 9.6e-3 * fbh_**-0.65 * np.exp(0.03 * (np.log(fbh_))**2)
        S2 = min(1, S2)
    else:
        S2 = 1

    S = S1 * S2
    
    return S
    

def dR_dm1dm2_kernel(fbh, mc, sbh, m1, m2, t, m2_ave, m_ave, mf_model, Use_S2, S1_method, Use_interp, S_Table):
    '''
    Differential comoving merger rate, unit : 1 / ( Gpc^3 * yr * msun^2)
    '''

    t0 = 4.3523939036782445e+17
    M = m1 + m2
    eta = m1*m2/M**2
    if Use_interp:
        S = np.interp(x = np.log(M), xp = S_Table['LnM'], fp = S_Table['S'], left = S_Table['S'][0], right = S_Table['S'][-1])
    else:
        S = S_Factor(fbh = fbh, sbh = sbh, mc = mc, M = M, m_ave = m_ave, m2_ave = m2_ave, mf_model = mf_model, t = t, 
                     Use_S2 = Use_S2, S1_method = S1_method)
    
    Psi_1 = Psi(mc = mc, sbh = sbh, m = m1, mf_model = mf_model)
    Psi_2 = Psi(mc = mc, sbh = sbh, m = m2, mf_model = mf_model)

    r1 = 1.6e6 * fbh**(53/37) * (t/t0)**(-34/37) * M**(-32/37) * eta**(-34/37)
    r2 = S * Psi_1 * Psi_2 / (m_ave**2)

    result = r1*r2
    
    return result

def Merger_Rate_Lite(fbh, mc, z, Use_S2, S1_method):

    t0 = 4.3523939036782445e+17
    t = p21ct.z2t(z)
    R0 = 3.14061016E6 * mc**(-32/37) * fbh**(53/37) * (t/t0)**(-34/37)
    S = S_Factor(fbh = fbh, sbh = 0, mc = mc, M = 2*mc, m_ave = mc, m2_ave = mc**2, mf_model = 2, t = t, Use_S2 = Use_S2, S1_method = S1_method)
    R = R0*S
    #print(t)
    return R

def get_dR_dm1dm2_vec(fbh, mc, sbh, z, mf_model, sbh_width, nm, Use_S2, show_status, Use_interp, S1_method, S_Tab_Len):
    '''
    Get a 2d array of dR_dm1dm2
    '''

    def get_m_vec():
        u = np.log(mc)
        x1 = u - sbh_width * sbh
        x2 = u + sbh_width * sbh
        x_vec = np.linspace(x1, x2, nm)
        m_vec = np.exp(x_vec)
        return m_vec
    
    m_ave = m_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    m2_ave = m2_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    t = p21ct.z2t(z = z)
    m_vec = get_m_vec()

    def Get_S_Table():
        M1 = 2 * np.min(m_vec)
        M2 = 2 * np.max(m_vec)
        x1 = np.log(M1)
        x2 = np.log(M2)
        LnM = np.linspace(x1, x2, S_Tab_Len)
        M_vec = np.exp(LnM)
        S_Tab = np.zeros(S_Tab_Len)
        t1 = PyLab.TimeNow()
        for idx in np.arange(0, S_Tab_Len):
            if show_status:
                print('Get_S_Table status : ', idx/S_Tab_Len)
            S_Tab[idx] = S_Factor(fbh = fbh, sbh = sbh, mc = mc, M = M_vec[idx], m_ave = m_ave, m2_ave = m2_ave, 
                                  mf_model = mf_model, t = t, Use_S2 = Use_S2, S1_method = S1_method)
        t2 = PyLab.TimeNow()
        if show_status:
            print('S Tab ready, time used:', t2 - t1)
        r = {'LnM' : LnM, 'S' : S_Tab}
        return r
    
    if Use_interp:
        S_Tab = Get_S_Table()
    else:
        S_Tab = []

    dRdm1dm2 = np.zeros((nm, nm))

    for mid_1 in np.arange(0, nm):
        if show_status:
            print('get_dR_dm1dm2_vec status : ', mid_1/nm)

        for mid_2 in np.arange(0, nm):

            m1 = m_vec[mid_1]
            m2 = m_vec[mid_2]

            if mid_1 <= mid_2:
                # m1 and m2 are symetric so we only need to do half of calculations
                dRdm1dm2[mid_1, mid_2] = dR_dm1dm2_kernel(fbh = fbh, mc = mc, sbh = sbh, m1 = m1, m2 = m2, t = t, m2_ave = m2_ave, m_ave = m_ave, mf_model = mf_model, 
                                                          Use_S2 = Use_S2, S1_method = S1_method, Use_interp = Use_interp, S_Table = S_Tab)
    # Now do the other half(ish)
    for mid_1 in np.arange(0, nm):
        for mid_2 in np.arange(0, nm):

            m1 = m_vec[mid_1]
            m2 = m_vec[mid_2]

            if mid_1 > mid_2:
                dRdm1dm2[mid_1, mid_2] = dRdm1dm2[mid_2, mid_1]
    
    result = {'m_vec' : m_vec, 'dRdm1dm2': dRdm1dm2}

    return result

def Merger_Rate_Scale_Factor(z, fbh, Use_S2):
    '''
    '''
    
    PyLab.Is_Scalar(x = z, ErrorMethod = 2)
    
    def S2T(z_):
        t0 = 4.3523939036782445e+17
        t = p21ct.z2t(z = z_)
        fbh_ = fbh*(t/t0)**0.44
        S2= 9.6e-3 * fbh_**-0.65 * np.exp(0.03 * (np.log(fbh_))**2)
        S2 = min(1.0, S2)
        if not Use_S2:
            S2 = 1.0
        tf = (t/t0)**(-34/37)
        r = tf * S2
        
        return r
    
    f0 = S2T(0)
    f = S2T(z)
    r = f/f0

    return r

def Get_dEGW_dv(v, eta, M):
    
    c3 = 299792458.0**3.0
    G = 6.6740831313131E-11
    msun = 1.98847E30

    a = [2.9740e-1, 5.9411e-1, 8.4845e-1, 5.0801e-1]
    b = [4.4810e-2, 8.9794e-2, 1.2848e-1, 7.7515e-2]
    c = [9.5560e-2, 1.9111e-1, 2.7299e-1, 2.2369e-2]
    
    v1 = a[0] * eta**2 + b[0] * eta + c[0]
    v2 = a[1] * eta**2 + b[1] * eta + c[1]
    v3 = a[2] * eta**2 + b[2] * eta + c[2]
    v4 = a[3] * eta**2 + b[3] * eta + c[3]
    
    v1 = v1 * c3/(np.pi * G * M * msun)
    v2 = v2 * c3/(np.pi * G * M * msun)
    v3 = v3 * c3/(np.pi * G * M * msun)
    v4 = v4 * c3/(np.pi * G * M * msun)
    
    if  v < v1:
        F = v**(-1/3)
    elif v1 <=v and v < v2:
        F = v**(2/3)/v1
    elif v2 <= v and v < v3:
        # There was a typo in 1707.01480, v should be v^2
        f1 = v**2 * v4**4
        f2 = v1 * v2**(4/3) * (4*(v - v2)**2 + v4**2)**2    
        F = f1/f2
    else:
        F = 0
    r = (np.pi * G)**(2/3) * (M*msun)**(5/3) * eta * F
    return r

def Get_OmGW_Lite(fbh, mc, v_vec, Use_S2, S1_method, nz, Precision, ncpu, show_status):
    
    if nz > 200:
        warnings.warn('nz is too high, this will be very slow')

    PyLab.Is_Scalar(v_vec, 1)
    
    def Find_dOmGW_dx(zp, v):
        RhoCrC2 = 7.730937688449169e-10
        Gpc3yr = 9.274526196872258e+83
        n = len(zp)
        dOmGW_dz = np.zeros(n)

        for idx in np.arange(0, n):
            z = zp[idx] - 1
            R = Merger_Rate_Lite(fbh = fbh, mc = mc, z = z, Use_S2 = Use_S2, S1_method = S1_method)
            vr = v*(1+z)
            dEGW_dvr = Get_dEGW_dv(v = vr, eta = 1/4, M = 2 * mc)/(1+z)
            # dEGW_dvr = Get_dEGW_dv(v = vr, eta = 1/4, M = 2 * mc)
            
            H = PyLab.Hubble(z=z)

            # dOmGW_dz[idx] = v/RhoCrC2 * R * dEGW_dvr / (H*(1+z))
            # print('debugging, a factor of (1+z) is missing')
            dOmGW_dz[idx] = v/RhoCrC2 * R * dEGW_dvr / (H*(1+z))
        
        dOmGW_dlnzp = zp*dOmGW_dz/Gpc3yr
        
        return dOmGW_dlnzp

    def Find_OmGW(v):
        '''
        Get OmGW for each v
        '''
        kernel = lambda zp: Find_dOmGW_dx(zp = zp, v = v)
        zp, fx = PyLab.Map(
            F = kernel,
            Start = 1,
            Width = 1,
            MinX = 0,
            MaxX = np.inf,
            nx = nz,
            Precision = Precision,
            Max_Iteration = 50,
            Use_log_x = 1,
            flat_count = 10)
        
        x = np.log(zp)
        r = np.trapz(x = x, y = fx)
        
        return r
    
    if ncpu < 2:
        nv = len(v_vec)
        r = np.zeros(nv)
        
        for idx in np.arange(0, nv):
            if show_status:
                print('Get_OmGW_Lite status:', idx/nv)
            r[idx] = Find_OmGW(v = v_vec[idx])
    else:
        r = Parallel(n_jobs=ncpu)(delayed(Find_OmGW)(v) for v in v_vec)
    
    return r

def dOmGW_dlnv_kernel_fast(v, nz, dR0_dm1dm2_info, z_info):
    ''' 
    '''

    # Converting to SI unit
    # RhoCrC2 = 8.6018282524e-27 * 299792458.0**2
    # yr = 365.25 * 24 * 3600
    # Gpc = 3.086E25
    # Gpc3yr = Gpc**3*yr

    RhoCrC2 = 7.730937688449169e-10
    Gpc3yr = 9.274526196872258e+83
    
    z_vec = z_info['z_vec']
    zp_vec = z_vec + 1.0
    
    m_vec = dR0_dm1dm2_info['m_vec']
    dR0dm1dm2_vec = dR0_dm1dm2_info['dRdm1dm2']
    nm = len(m_vec)
    dOmGW_dlnm1_dlnm2_dz = np.zeros(nm)
    dOmGW_dlnm1_dz = np.zeros(nm)
    dOmGW_dz = np.zeros(nz)

    for zid in np.arange(0, nz):
        # print(zid/nz)
        z_ = z_vec[zid]
        H = z_info['H_vec'][zid]
        z_factor = z_info['scale_factor'][zid]
        vr = v * (1+z_)

        for mid_1 in np.arange(0, nm):
            for mid_2 in np.arange(0, nm):

                m1 = m_vec[mid_1]
                m2 = m_vec[mid_2]
                dR0dm1dm2 = dR0dm1dm2_vec[mid_1, mid_2]
                dRdm1dm2 = dR0dm1dm2 * z_factor
                dR_dlnm1_dlnm2 = dRdm1dm2 * m1 * m2
                M = m1 + m2
                eta = m1*m2/M**2
                dEGW_dvr = Get_dEGW_dv(v = vr, eta = eta, M = M)/(1+z_)
                # dEGW_dvr = Get_dEGW_dv(v = vr, eta = eta, M = M)
                
                dOmGW_dlnm1_dlnm2_dz[mid_2] = v/RhoCrC2 * dR_dlnm1_dlnm2 * dEGW_dvr/((1+z_)*H)
            
            
            dOmGW_dlnm1_dz[mid_1] = np.trapz(x = np.log(m_vec), y = dOmGW_dlnm1_dlnm2_dz)
        dOmGW_dz[zid] = np.trapz(x = np.log(m_vec), y = dOmGW_dlnm1_dz)
    
    dOmGW_dlnzp = zp_vec * dOmGW_dz
    dOmGW_dlnv = np.trapz(x = np.log(zp_vec), y = dOmGW_dlnzp)

    dOmGW_dlnv = dOmGW_dlnv/Gpc3yr
    
    return dOmGW_dlnv

def Get_dOmGW_dlnzp(zp, v, fbh, Use_S2, Table):

    dR0_dm1dm2_vec = Table['dR0_dm1dm2_vec']
    m_vec = Table['m_vec']
    z = zp-1
    
    nz = len(z)
    nm = len(m_vec)
    lnm = np.log(m_vec)
    fz = np.zeros(nz)
    f1 = np.zeros(nm)
    f2 = np.zeros(nm)
    
    for idx in np.arange(0, nz):
        z_ = z[idx]
        Scale_Factor = Merger_Rate_Scale_Factor(z = z_, fbh = fbh, Use_S2 = Use_S2)
        vr = v*(1+z_)

        for mid_1 in np.arange(0, nm):
            for mid_2 in np.arange(0, nm):

                m1 = m_vec[mid_1]
                m2 = m_vec[mid_2]
                dR0_dm1dm2 = dR0_dm1dm2_vec[mid_1, mid_2]
                dR_dlnm1dlnm2 = Scale_Factor * dR0_dm1dm2 * m1 * m2

                M = m1 + m2
                eta = m1*m2/M**2
                # dEGW_dvr = Get_dEGW_dv(v = vr, eta = eta, M = M)
                dEGW_dvr = Get_dEGW_dv(v = vr, eta = eta, M = M)/(1+z_)

                f2[mid_2] = dR_dlnm1dlnm2 * dEGW_dvr
            
            f1[mid_1] = np.trapz(x = lnm, y = f2)
        fz[idx] = np.trapz(x = lnm, y = f1)
    
    RhoCrC2 = 7.730937688449169e-10
    Gpc3yr = 9.274526196872258e+83
    H = PyLab.Hubble(z = z)
    
    # OmGW = \int dz Fz 
    Fz = v/RhoCrC2 * fz / (H * (1+z))
    Fz = Fz/Gpc3yr
    
    # but we want to do this in log axis,
    # OmGW = \int dlnzp Fx
    Fx = zp*Fz
    
    return Fx


def Get_OmGW_Pro(v_vec, nz, fbh, Use_S2, mc, sbh, mf_model, sbh_width, nm, show_status, Use_interp, S1_method, S_Tab_Len, Precision, ncpu):
    '''
    Find OmGW for each individual v
    '''
    
    Tables = get_dR_dm1dm2_vec(
        fbh = fbh, 
        mc = mc, 
        sbh = sbh, 
        z = 0, 
        mf_model = mf_model, 
        sbh_width = sbh_width, 
        nm = nm, 
        Use_S2 = Use_S2, 
        show_status = show_status,
        Use_interp = Use_interp, 
        S1_method = S1_method, 
        S_Tab_Len = S_Tab_Len)
    
    Table = {'m_vec' : Tables['m_vec'], 'dR0_dm1dm2_vec' : Tables['dRdm1dm2']}
    
    # check Psi convergence
    m_vec = Tables['m_vec']
    mf = Psi(mc = mc, sbh = sbh, m = m_vec, mf_model = mf_model)
    Norm = np.trapz(x = np.log(m_vec), y = mf)
    if abs(Norm - 1) > 0.1:
        raise Exception('Mass function is not convergent, please increase precision by setting sbh_width and nm')
    
    def kernel(v):
        # Get OmGW for each v
        fun = lambda zp: Get_dOmGW_dlnzp(zp = zp, v = v, fbh = fbh, Use_S2 = Use_S2, Table = Table)
        
        zp, fx = PyLab.Map(
            F = fun,
            Start = 1,
            Width = 1,
            MinX = 0,
            MaxX = np.inf,
            nx = nz,
            Precision = Precision,
            Max_Iteration = 50,
            Use_log_x = 1,
            Print_debug_MSG = 0,
            flat_count = 10)
        x = np.log(zp)
        r = np.trapz(x = x, y = fx)
        return r
    
    if ncpu < 2:
        nv = len(v_vec)
        r = np.zeros(nv)
        for idx in np.arange(0, nv):
            if show_status:
                print('Get_OmGW_Pro status: ', idx/nv)
            r[idx] = kernel(v_vec[idx])
    else:
        r = Parallel(n_jobs = ncpu)(delayed(kernel)(v) for v in v_vec)
    
    return r

def Merger_Rate(
        fbh = 1e-2, 
        mc = 20, 
        sbh = 0.6, 
        z = 0, 
        mf_model = 1, 
        sbh_width = 10, 
        nm = 200,
        Use_S2 = 0, 
        S1_method = 0,
        Use_interp = 1,
        S_Tab_Len = 200,
        show_status = 0):
    
    if mf_model == 2:
        if PyLab.Is_Scalar(z):
            r = Merger_Rate_Lite(fbh = fbh, mc = mc, z = z, Use_S2 = Use_S2, S1_method = S1_method)
        else:
            nz = len(z)
            r = np.zeros(nz)
            for idx in np.arange(0, nz):
                r[idx] = Merger_Rate_Lite(fbh = fbh, mc = mc, z = z[idx], Use_S2 = Use_S2, S1_method = S1_method)
        return r
    
    # First get z=0 result for scaling
    tab = get_dR_dm1dm2_vec(fbh = fbh, mc = mc, sbh = sbh, z = 0, mf_model = mf_model, sbh_width = sbh_width, nm = nm, Use_S2 = Use_S2, 
                            show_status = show_status, Use_interp = Use_interp, S1_method = S1_method, S_Tab_Len = S_Tab_Len)

    m = tab['m_vec']
    dRdm1dm2 = tab['dRdm1dm2']
    x = np.log(m)
    dR_dlnm1_dlnm2 = np.zeros(nm)
    dR_dlnm1 = np.zeros(nm)

    for idx_1 in np.arange(0, nm):
        for idx_2 in np.arange(0, nm):

            m1 = m[idx_1]
            m2 = m[idx_2]
            
            dR_dlnm1_dlnm2[idx_2] = dRdm1dm2[idx_1, idx_2] * m1 * m2
            
        dR_dlnm1[idx_1] = np.trapz(x = x, y = dR_dlnm1_dlnm2)
    
    R0 = np.trapz(x = x, y = dR_dlnm1) # for z=0
    
    if PyLab.Is_Scalar(z):
        r = R0 * Merger_Rate_Scale_Factor(z = z, fbh = fbh, Use_S2 = Use_S2)
    else:
        nz = len(z)
        r = np.zeros(nz)
        for idx in np.arange(0, nz):
            r[idx] = R0 * Merger_Rate_Scale_Factor(z = z[idx], fbh = fbh, Use_S2 = Use_S2)
        
    return r

def Get_dOmGW_dlnv(
        fbh = 1e-2, 
        mc = 30,
        sbh = 1.0, 
        v = np.logspace(-6, 4, 50), 
        mf_model = 0, 
        sbh_width = 10, 
        nm = 200,
        nz = 100,
        zmax = 50,
        show_status = 0,
        Use_interp = 1,
        S1_method = 0,
        Fast = 0,
        S_Tab_Len = 200,
        Use_S2 = 1,
        Precision = 1e-2,
        ncpu = 1):
    '''
    What to feed to dOmGW_dlnv_kernel_fast:
    z_vec
    H
    z_factor
    do Gpcyr outside
    '''
    t1 = PyLab.TimeNow()
    if mf_model == 2:
        r = Get_OmGW_Lite(fbh = fbh, mc = mc, v_vec = v, Use_S2 = Use_S2, 
                          S1_method = S1_method, nz = nz,Precision = Precision, ncpu = ncpu, show_status = show_status)
        t2 = PyLab.TimeNow()
        if show_status:
            print('Get_dOmGW_dlnv: all done, time used:', t2 - t1)
        return r
    
    if not Fast:
        r = Get_OmGW_Pro(
            fbh = fbh, 
            mc = mc, 
            sbh = sbh, 
            v_vec = v, 
            mf_model = mf_model, 
            nz = nz, 
            Use_S2 = Use_S2, 
            sbh_width = sbh_width, 
            nm = nm, 
            show_status = show_status, 
            Use_interp = Use_interp, 
            S1_method = S1_method, 
            S_Tab_Len = S_Tab_Len, 
            Precision = Precision, 
            ncpu = ncpu)
        t2 = PyLab.TimeNow()
        if show_status:
            print('Get_dOmGW_dlnv: all done, time used:', t2 - t1)
        return r
    
    # ---- fast mode ----
    def get_z_info():
        zp_vec = np.logspace(0.0, np.log10(1+zmax), nz)
        z_vec = zp_vec - 1
        H_vec = PyLab.Hubble(z_vec)
        scale_factor = np.zeros(nz)
        for idx in np.arange(0, nz):
            scale_factor[idx] = Merger_Rate_Scale_Factor(z = z_vec[idx], fbh = fbh, Use_S2 = Use_S2)
        r = {'z_vec' : z_vec, 'H_vec' : H_vec, 'scale_factor': scale_factor}
        return r
    
    nv = len(v)
    r = np.zeros(nv)
    z_info = get_z_info()
    
    if show_status:
        print('Computing merger rate template')
    
    dR0_dm1dm2_info = get_dR_dm1dm2_vec(fbh = fbh, mc = mc, sbh = sbh, z = 0, mf_model = mf_model, sbh_width = sbh_width, nm = nm, Use_S2 = Use_S2, 
                                        show_status = show_status, Use_interp = Use_interp, S1_method = S1_method, S_Tab_Len = S_Tab_Len)
    
    # At the very least, Psi should be convergent

    m_vec = dR0_dm1dm2_info['m_vec']
    mf = Psi(mc = mc, sbh = sbh, m = m_vec, mf_model = mf_model)
    Norm = np.trapz(x = np.log(m_vec), y = mf)
    if abs(Norm - 1) > 0.1:
        raise Exception('Mass function is not convergent, please increase precision by setting sbh_width and nm')
    
    def model(x):
        r = dOmGW_dlnv_kernel_fast(v = x, nz = nz, dR0_dm1dm2_info = dR0_dm1dm2_info, z_info = z_info)
        PyLab.SaySomething()
        return r
    
    if ncpu == 1:
        for idx in np.arange(0, nv):
            if show_status:
                print('dOmGW_dlnv status:', idx/nv)
            v_ = v[idx]
            if v_<0:
                r[idx] = np.nan
            else:
                r[idx] = model(v_)
    else:
        r = Parallel(n_jobs = ncpu)(delayed(model)(x) for x in v)
    
    if (True in r<0) or (True in np.isnan(r)):
        print('fbh = ', fbh, ', mc = ', mc, ', sbh = ', sbh)
        raise Exception('Something went wrong, OmGW is negative')
    t2 = PyLab.TimeNow()
    if show_status:
        print('Get_dOmGW_dlnv: all done, time used:', t2 - t1)
    return r
