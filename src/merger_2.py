import numpy as np
import scipy.special as sp
import p21c_tools as p21ct
import mpmath
import PyLab
import scipy.io
from joblib import Parallel, delayed

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

clear
fbh_vec = logspace(-20, 0, 10000);
sigma_M = 0.005939410849858249;
for idx=1:length(fbh_vec)
    U_vec(idx) = kummerU(21/74, 1/2, 5*fbh_vec(idx)^2/(6*sigma_M^2));
end
save /Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template_2.mat fbh_vec U_vec sigma_M
'''

kummerU_template = scipy.io.loadmat('/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template.mat')
kummerU_template = scipy.io.loadmat('/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template_2.mat')
sigma_M = 0.004
#sigma_M = 0.005939410849858249
S1_nm = 300
Use_Fast_S1 = 0
S1_nv = 200
sbh_width_s1 = 6
t0 = 4.3523939036782445e+17
s_table_size = 500

def Psi(mc, sbh, m, mf_model):
    '''
    dfbh/dm/fbh
    ----inputs----
    mf_model : MF model
        0 : log-normal in rho
        1 : log-normal in nbh
        2 : monochromatic
    '''
    x = np.log(m)
    u = np.log(mc)
    chi2 = ((x-u)/sbh)**2
    result = np.exp(-chi2/2)/(sbh*m*np.sqrt(2 * np.pi))
    return result

def Get_P(mc, sbh, mf_model):
    n = 1000
    sbh_width = 6
    u = np.log(mc)
    x1 = u - sbh * sbh_width
    x2 = u + sbh * sbh_width
    x_vec = np.linspace(x1, x2, n)
    m_vec = np.exp(x_vec)
    mf = Psi(mc = mc, sbh = sbh, m = m_vec, mf_model = mf_model)
    P_inv = np.trapz(x = x_vec, y = mf)
    P = 1/P_inv
    return P, mf, m_vec

def m_ave_fun(mc, sbh, mf_model):
    '''
    Mean of m
    '''
    # result = mc*np.exp(sbh**2/2)
    
    P, mf, m = Get_P(mc = mc, sbh = sbh, mf_model = mf_model)
    x = np.log(m)
    fun = P * m * mf
    r = np.trapz(x = x, y = fun)
    return r

def m2_ave_fun(mc, sbh, mf_model):
    '''
    Mean of m^2
    '''
    result = mc**2 * np.exp(2 * sbh**2)
    return result

def Get_S1_Fast(fbh, sbh, M, m_ave, m2_ave):
    '''
    '''

    sigma_M = kummerU_template['sigma_M'][0][0]
    Ny = M*fbh/(m_ave*(fbh + sigma_M))
    
    # Get C
    if sbh > 2:
        raise Exception('Approximation for C will be compromised for sbh > 2!')
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

def Get_S1(fbh, mc, sbh, M, mf_model):
    
    m_ave = m_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)

    if Use_Fast_S1:
        m2_ave = m2_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
        S1 = Get_S1_Fast(fbh = fbh, sbh = sbh, M = M, m_ave = m_ave, m2_ave = m2_ave)
        return S1

    Ny = M*fbh/(m_ave * (fbh + sigma_M))
    
    def int_m_kernel(v):
        u = np.log(mc)
        x1 = u - sbh_width_s1 * sbh
        x2 = u + sbh_width_s1 * sbh
        x_vec = np.linspace(x1, x2, S1_nm)
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
        
        p3 = Ny*m_ave*int_m_vec(v)
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
        nx = S1_nv,
        Precision = 5e-3,
        Max_Iteration = 50,
        Use_log_x = 0,
        flat_count = 20)
    
    r = np.trapz(x = np.log(v), y = v * fv)
    
    return r

def S_factor(fbh, mc, sbh, M, m_ave, mf_model, t, Type):
    '''
    Type : 
        0 - S1
        1 - S1*S2
        2 - S2
    '''
    
    if Type in [1, 2]:
        fbh_ = fbh*(t/t0)**0.44
        S2= 9.6e-3 * fbh_**-0.65 * np.exp(0.03 * (np.log(fbh_))**2)
        S2 = min(1, S2)
    if Type in [0, 1]:
        S1 = Get_S1(fbh = fbh, mc = mc, sbh = sbh, M = M, mf_model = mf_model)
    
    if Type == 0:
        r = S1
    elif Type == 1:
        r = S1*S2
    else:
        r = S2
    
    return r
    
def dR_dm1dm2_z0_kernel(fbh, mc, sbh, m1, m2, m_ave, use_interp, mf_model, Use_S2, lnM_axis, S_axis):
    '''
    Differential comoving merger rate at z = 0, unit : 1 / ( Gpc^3 * yr * msun^2)
    '''

    M = m1 + m2
    eta = m1*m2/M**2
    psi_1 = Psi(mc = mc, sbh = sbh, m = m1, mf_model = mf_model)
    psi_2 = Psi(mc = mc, sbh = sbh, m = m2, mf_model = mf_model)

    dR0_dm1_dm2 = 1.6e6 * fbh**(53/37) * eta**(-34/37) * M**(-32/37) * psi_1 * psi_2
    
    if use_interp == 0:
        S = S_factor(fbh = fbh, mc = mc, sbh = sbh, M = M, m_ave = m_ave, mf_model = mf_model, t = t0, Type = Use_S2)
    else:
        S = np.interp(x = np.log(M), xp = lnM_axis, fp = S_axis)
    r = dR0_dm1_dm2 * S
    return r    

def Scale_Factor(z, fbh, Use_S2):
    t = p21ct.z2t(z = z)
    t_factor = (t/t0)**(-34/37)
    if Use_S2:
        S2_0 = S_factor(fbh = fbh, t = t0, Type = 2, mc = [], sbh = [], M = [], m_ave = [], mf_model = [])
        S2 = S_factor(fbh = fbh, t = t, Type = 2, mc = [], sbh = [], M = [], m_ave = [], mf_model = [])
        t_factor = t_factor * S2/S2_0
    
    return t_factor

def get_dR_dm1dm2_vec(fbh, mc, sbh, z, sbh_width, nm, use_interp, show_status, Use_S2, mf_model):
    '''
    Get a 2d array of dR_dm1dm2
    '''
    
    t_factor = Scale_Factor(z = z, fbh = fbh, Use_S2 = Use_S2)
    t = p21ct.z2t(z = z)
    m_ave = m_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    
    def get_m_vec():
        u = np.log(mc)
        x1 = u - sbh_width * sbh
        x2 = u + sbh_width * sbh
        x_vec = np.linspace(x1, x2, nm)
        m_vec = np.exp(x_vec)
        return m_vec
    
    def get_S_table(m_vec):
        # Get an interpolation table for S
        n = min(s_table_size, nm**2) # should be enough
        lnm1 = np.log(2 * np.min(m_vec))
        lnm2 = np.log(2 * np.max(m_vec))
        lnM_axis = np.linspace(lnm1, lnm2, n)
        M_vec = np.exp(lnM_axis)
        S_axis = np.zeros(n)
        
        for idx in np.arange(0, n):
            if show_status:
                print('get_S_table status : ', idx/n)
            S_axis[idx] = S_factor(fbh = fbh, mc = mc, sbh = sbh, M = M_vec[idx], m_ave = m_ave, mf_model = mf_model, t = t, Type = Use_S2)
        return lnM_axis, S_axis

    m_vec = get_m_vec()
    if use_interp:
        lnM_axis, S_axis = get_S_table(m_vec)
    else:
        lnM_axis, S_axis = [], []

    dRdm1dm2 = np.zeros((nm, nm))
    # m1 and m2 are symetric so we only need to do half of calculations
    for mid_1 in np.arange(0, nm):
        if show_status:
            print('get_dR_dm1dm2_vec status : ', mid_1/nm)
        for mid_2 in np.arange(0, nm):

            m1 = m_vec[mid_1]
            m2 = m_vec[mid_2]
            if mid_1 <= mid_2:
                dRdm1dm2_z0 = dR_dm1dm2_z0_kernel(fbh = fbh, mc = mc, sbh = sbh, m1 = m1, m2 = m2, m_ave = m_ave, use_interp = use_interp,
                                                  mf_model = mf_model, Use_S2 = Use_S2, lnM_axis = lnM_axis, S_axis = S_axis)

                dRdm1dm2[mid_1, mid_2] = t_factor * dRdm1dm2_z0
    
    for mid_1 in np.arange(0, nm):
        for mid_2 in np.arange(0, nm):
            if mid_1 > mid_2:
                dRdm1dm2[mid_1, mid_2] = dRdm1dm2[mid_2, mid_1]
    
    result = {'m_vec' : m_vec, 'dRdm1dm2': dRdm1dm2}

    return result

def Merger_Rate(
        fbh = 1e-2, 
        mc = 20, 
        sbh = 0.6, 
        z = 0, 
        sbh_width = 10, 
        nm = 200, 
        use_interp = 1, 
        Use_S2 = 0,
        mf_model = 0,
        show_status = 0):
    
    # First get z=0 result for scaling
    tab = get_dR_dm1dm2_vec(fbh = fbh, mc = mc, sbh = sbh, z = 0, sbh_width = sbh_width, 
                            nm = nm, use_interp = use_interp, show_status = show_status, Use_S2 = Use_S2, mf_model = mf_model)

    m = tab['m_vec']
    dRdm1dm2 = tab['dRdm1dm2']
    x = np.log(m)
    dR_dlnm1 = np.zeros(nm)
    dR_dlnm1_dlnm2 = np.zeros(nm)

    for idx_1 in np.arange(0, nm):
        for idx_2 in np.arange(0, nm):

            m1 = m[idx_1]
            m2 = m[idx_2]
            
            dR_dlnm1_dlnm2[idx_2] = dRdm1dm2[idx_1, idx_2] * m1 * m2
            
        dR_dlnm1[idx_1] = np.trapz(x = x, y = dR_dlnm1_dlnm2)
    
    R0 = np.trapz(x = x, y = dR_dlnm1)
    
    if PyLab.Is_Scalar(z):
        t_factor = Scale_Factor(z = z, fbh = fbh, Use_S2 = Use_S2)
        r = R0 * t_factor
    else:
        nz = len(z)
        r = np.zeros(nz)
        for idx in np.arange(0, nz):
            t_factor = Scale_Factor(z = z[idx], fbh = fbh, Use_S2 = Use_S2)
            r[idx] = R0 * t_factor
        
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
    elif v2 <=v and v < v3:
        # There was a typo in 1707.01480, v should be v^2
        f1 = v**2 * v4**4
        f2 = v1 * v2**(4/3) * (4*(v - v2)**2 + v4**2)**2
        F = f1/f2
    else:
        F = 0
    
    
    r = (np.pi * G)**(2/3) * (M*msun)**(5/3) * eta * F
    return r


def dOmGW_dlnv_kernel(v, nz, R0, z_info):
    ''' 
    actually even z dimension can also be simplified from outside,
    but this is supposed to be super fast anyway
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
    
    m_vec = R0['m_vec']
    dR0dm1dm2_vec = R0['dRdm1dm2']
    nm = len(m_vec)
    dOmGW_dlnm1_dlnm2_dz = np.zeros(nm)
    dOmGW_dlnm1_dz = np.zeros(nm)
    dOmGW_dz = np.zeros(nz)

    for zid in np.arange(0, nz):
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
                dEGW_dv = Get_dEGW_dv(v = vr, eta = eta, M = M)
                dEGW_dvr = dEGW_dv/(1+z_)
                
                # dOmGW_dlnm1_dlnm2_dz[mid_2] = v/RhoCrC2 * dR_dlnm1_dlnm2 * dEGW_dv/((1+z_)*H)
                dOmGW_dlnm1_dlnm2_dz[mid_2] = v/RhoCrC2 * dR_dlnm1_dlnm2 * dEGW_dvr/((1+z_)*H)
            
            
            dOmGW_dlnm1_dz[mid_1] = np.trapz(x = np.log(m_vec), y = dOmGW_dlnm1_dlnm2_dz)
        dOmGW_dz[zid] = np.trapz(x = np.log(m_vec), y = dOmGW_dlnm1_dz)
    
    dOmGW_dlnzp = zp_vec * dOmGW_dz
    dOmGW_dlnv = np.trapz(x = np.log(zp_vec), y = dOmGW_dlnzp)

    dOmGW_dlnv = dOmGW_dlnv/Gpc3yr
    
    return dOmGW_dlnv