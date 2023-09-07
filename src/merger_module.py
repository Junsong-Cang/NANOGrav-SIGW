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
save /Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template.mat fbh_vec U_vec
'''

kummerU_template = scipy.io.loadmat('/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src/kummerU_template.mat')

def Psi(mc, sbh, m, mf_model):
    '''
    dnbh/dlnm/nbh
    ----inputs----
    mf_model : MF model
        0 : log-normal in rho
        1 : log-normal in nbh
        2 : log10-normal in rho
        3 : monochromatic
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

def m_ave_fun(mc, sbh, mf_model):
    '''
    Mean of m
    '''
    sbh_width = 10
    nm = 10000

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
    nm = 10000

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

def dR_dm1dm2_kernel(fbh, mc, sbh, m1, m2, t, m2_ave, m_ave, mf_model):
    '''
    Differential comoving merger rate, unit : 1 / ( Gpc^3 * yr * msun^2)
    I want to avoid repeated computation at the cost of simplicity, m_ave and m2_ave can be computed by:
    m_ave = m_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    m2_ave = m2_ave_fun(mc = mc, sbh = sbh, mf_model = mf_model)
    '''
    debug = 0
    t0 = 4.3523939036782445e+17
    M = m1 + m2
    eta = m1*m2/M**2
    
    def Get_S():
        # do this in seperate block to avoid name space pollution and repeated calculation

        # 1: N(y)
        sigma_M = 0.004
        Ny = M*fbh/(m_ave*(fbh + sigma_M))
        
        # 2: C
        if sbh >2:
            raise Exception('Approximation for c will be compromised!')
        c1 = (fbh/(sigma_M*m_ave))**2 * m2_ave
        c2a = sp.gamma(29/37)/np.sqrt(np.pi)
        # I donno why but scipy hyp1f1 is very strange!
        if debug:
            U = sp.hyp1f1(21/74, 1/2, 5*fbh**2/(6*sigma_M**2))
        else:
            log_fbh_axis = np.log(kummerU_template['fbh_vec'][0,:])
            U_axis = kummerU_template['U_vec'][0,:]
            U = np.interp(x = np.log(fbh), xp = log_fbh_axis, fp = U_axis)
        
        c2 = (c2a*U)**(-74/21)-1
        C = c1/c2

        # 3 : S1
        s1a = 1.42*np.exp(-Ny)
        s1b = m2_ave/(m_ave**2)/(Ny + C)
        s1c = (sigma_M/fbh)**2
        S1 = s1a * (s1b + s1c)**(-21/74)
        if np.isnan(S1):
            #print(s1a, s1b , s1c, (s1b + s1c)**(-21/74), C, U)
            print(U, 5*fbh**2/(6*sigma_M**2))
            
            raise Exception(' found -----')


        # 4 : s2
        fbh_ = fbh*(t/t0)**0.44
        s2a = 1.0
        s2b= 9.6e-3 * fbh_**-0.65 * np.exp(0.03 * (np.log(fbh_))**2)
        S2 = min(s2a, s2b)

        # finally get S
        S = S1 * S2
        return S
    
    S = Get_S()
    
    p1 = Psi(mc = mc, sbh = sbh, m = m1, mf_model = mf_model)
    p2 = Psi(mc = mc, sbh = sbh, m = m2, mf_model = mf_model)

    r1 = 1.6e6 * fbh**(53/37) * (t/t0)**(-34/37) * M**(-32/37) * eta**(-34/37)
    r2 = S*p1*p2/(m_ave**2)

    result = r1*r2
    
    return result

def get_dR_dm1dm2_vec(fbh, mc, sbh, z, mf_model, sbh_width, nm, show_status):
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
    t = p21ct.cosmic_age(z = z)
    m_vec = get_m_vec()
    
    dRdm1dm2 = np.zeros((nm, nm))

    for mid_1 in np.arange(0, nm):
        if show_status:
            print('get_dR_dm1dm2_vec status : ', mid_1/nm)

        for mid_2 in np.arange(0, nm):

            m1 = m_vec[mid_1]
            m2 = m_vec[mid_2]
            dRdm1dm2[mid_1, mid_2] = dR_dm1dm2_kernel(fbh = fbh, mc = mc, sbh = sbh, m1 = m1, m2 = m2, 
                                                      t = t, mf_model = mf_model, m_ave = m_ave, m2_ave = m2_ave)
    
    result = {'m_vec' : m_vec, 'dRdm1dm2': dRdm1dm2}

    return result

def Merger_Rate_Scale_Factor(z = 10.1, fbh = 1e-2):
    '''
    
    '''
    
    PyLab.Is_Scalar(x = z, ErrorMethod = 2)
    
    def S2T(z_):
        t0 = 4.3523939036782445e+17
        t = p21ct.cosmic_age(z = z_)
        fbh_ = fbh*(t/t0)**0.44
        s2a = 1.0
        s2b= 9.6e-3 * fbh_**-0.65 * np.exp(0.03 * (np.log(fbh_))**2)
        S2 = min(s2a, s2b)

        tf = (t/t0)**(-34/37)
        r = tf * S2
        
        return r
    
    f0 = S2T(0)
    f = S2T(z)
    r = f/f0

    return r

def Merger_Rate(fbh = 1e-2, mc = 20, sbh = 0.6, z = 0, mf_model = 1, sbh_width = 10, nm = 1000, show_status = 0):
    
    # First get z=0 result for scaling

    tab = get_dR_dm1dm2_vec(fbh = fbh, mc = mc, sbh = sbh, z = 0, mf_model = mf_model, sbh_width = sbh_width, nm = nm, show_status = show_status)

    m = tab['m_vec']
    dRdm1dm2 = tab['dRdm1dm2']
    x = np.log(m)
    f1 = np.zeros(nm)
    f2 = np.zeros(nm)

    for idx_1 in np.arange(0, nm):
        for idx_2 in np.arange(0, nm):

            m1 = m[idx_1]
            m2 = m[idx_2]
            
            f2[idx_2] = dRdm1dm2[idx_1, idx_2] * m1 * m2
            
        f1[idx_1] = np.trapz(x = x, y = f2)
    
    R0 = np.trapz(x = x, y = f1)
    
    if PyLab.Is_Scalar(z):
        r = R0 * Merger_Rate_Scale_Factor(z = z, fbh = fbh)
    else:
        nz = len(z)
        r = np.zeros(nz)
        for idx in np.arange(0, nz):
            r[idx] = R0 * Merger_Rate_Scale_Factor(z = z[idx], fbh = fbh)
        
    return r

def Get_dEGW_dv(v, eta, M):

    c3 = 299792458.0**3.0
    G = 6.6740831313131E-11
    msun = 1.98847E30

    a = [2.9740e-1, 5.9411e-1, 8.4845e-1, 5.0801e-1]
    b = [4.4810e-2, 8.9794e-2, 1.2848e-1, 7.7515e-2]
    c = [9.5560e-2, 1.9111e-1, 2.7299e-1, 2.2369e-2]
    
    v1 = a[0] * eta**2 + b[0] * eta * c[0]
    v2 = a[1] * eta**2 + b[1] * eta * c[1]
    v3 = a[2] * eta**2 + b[2] * eta * c[2]
    v4 = a[3] * eta**2 + b[3] * eta * c[3]

    v1 = v1 * c3/(np.pi * G * M * msun)
    v2 = v2 * c3/(np.pi * G * M * msun)
    v3 = v3 * c3/(np.pi * G * M * msun)
    v4 = v4 * c3/(np.pi * G * M * msun)

    if v < v1:
        F = v**(-1/3)
    elif v1 <=v and v < v2:
        F = v**(2/3)/v1
    elif v2 <=v and v < v3:
        # There was a typo in 1707.01480, v should vbe v^2
        f1 = v**2 * v4**4
        f2 = v1 * v2**(4/3) * (4*(v - v2)**2 + v4**2)**2
        F = f1/f2
    else:
        F = 0
    
    r = (np.pi * G)**(2/3) * (M*msun)**(5/3) * eta * F
    
    return r

def dOmGW_dlnv_kernel(v, nz, zmax, fbh, R0, z_info):
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
                # dEGW_dv = Get_dEGW_dv(v = vr, eta = eta, M = M)/(1+z_)
                dEGW_dv = Get_dEGW_dv(v = vr, eta = eta, M = M)
                dOmGW_dlnm1_dlnm2_dz[mid_2] = v/RhoCrC2 * dR_dlnm1_dlnm2 * dEGW_dv/((1+z_)*H)
            
            dOmGW_dlnm1_dz[mid_1] = np.trapz(x = np.log(m_vec), y = dOmGW_dlnm1_dlnm2_dz)
        dOmGW_dz[zid] = np.trapz(x = np.log(m_vec), y = dOmGW_dlnm1_dz)
    
    dOmGW_dlnzp = zp_vec * dOmGW_dz
    dOmGW_dlnv = np.trapz(x = np.log(zp_vec), y = dOmGW_dlnzp)

    dOmGW_dlnv = dOmGW_dlnv/Gpc3yr
    
    return dOmGW_dlnv

def Get_dOmGW_dlnv(
        fbh = 1e-2, 
        mc = 30,
        sbh = 1.0, 
        v = np.logspace(-6, 4, 50), 
        mf_model = 0, 
        sbh_width = 10, 
        nm = 200,
        nz = 200,
        zmax = 50,
        show_status = 0,
        ncpu = 1):
    '''
    WHat to feed to dOmGW_dlnv_kernel:
    z_vec
    H
    z_factor
    do Gpcyr outside
    '''
    def get_z_info():
        zp_vec = np.logspace(0.0, np.log10(1+zmax), nz)
        z_vec = zp_vec - 1
        H_vec = PyLab.Hubble(z_vec)
        scale_factor = np.zeros(nz)
        for idx in np.arange(0, nz):
            scale_factor[idx] = Merger_Rate_Scale_Factor(z = z_vec[idx], fbh = fbh)
        r = {'z_vec' : z_vec, 'H_vec' : H_vec, 'scale_factor': scale_factor}
        return r
    
    z_info = get_z_info()

    if show_status:
        print('Computing merger rate template')
    
    R0 = get_dR_dm1dm2_vec(fbh = fbh, mc = mc, sbh = sbh, z = 0, mf_model = mf_model, sbh_width = sbh_width, nm = nm, show_status = show_status)
    nv = len(v)
    r = np.zeros(nv)
    
    def model(x):
        r = dOmGW_dlnv_kernel(v = x, nz = nz, zmax = zmax, fbh = fbh, R0 = R0, z_info = z_info)
        return r
    
    if ncpu == 1:
        for idx in np.arange(0, nv):
            if show_status:
                print('dOmGW_dlnv status:', idx/nv)
            v_ = v[idx]
            r[idx] = model(v_)
    else:
            r = Parallel(n_jobs = ncpu)(delayed(model)(x) for x in v)
    return r

'''
def dR_dm1dm2_kernel_tmp(fbh = 0.05, mc = 30):
    debug = 0
    M = 2*mc
    m2_ave = mc**2
    m_ave = mc
    
    def Get_S():
        # do this in seperate block to avoid name space pollution and repeated calculation

        # 1: N(y)
        sigma_M=np.sqrt(0.005*1.137*1.137)
        Ny = M*fbh/(m_ave*(fbh + sigma_M))
        
        # 2: C
        c1 = (fbh/(sigma_M*m_ave))**2 * m2_ave
        c2a = sp.gamma(29/37)/np.sqrt(np.pi)
        # I donno why but scipy hyp1f1 is very strange!
        if debug:
            U = sp.hyp1f1(21/74, 1/2, 5*fbh**2/(6*sigma_M**2))
        else:
            log_fbh_axis = np.log(kummerU_template['fbh_vec'][0,:])
            U_axis = kummerU_template['U_vec'][0,:]
            U = np.interp(x = np.log(fbh), xp = log_fbh_axis, fp = U_axis)
        
        c2 = (c2a*U)**(-74/21)-1
        C = c1/c2

        # 3 : S1
        s1a = 1.42*np.exp(-Ny)
        s1b = m2_ave/(m_ave**2)/(Ny + C)
        s1c = (sigma_M/fbh)**2
        S1 = s1a * (s1b + s1c)**(-21/74)
        print(S1)

        return S1
    Get_S()

dR_dm1dm2_kernel_tmp()
'''
