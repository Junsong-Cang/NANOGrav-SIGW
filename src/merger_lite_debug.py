def Merger_Rate_Gouttenoire(fbh, m, z):
    
    msun = const['msun']
    RhoCr = cosmo['RhoCr']
    OmC = cosmo['OmC']
    
    mbh = m*msun
    RhoBH = OmC * fbh * RhoCr * (1+z)**3
    nbh = RhoBH / mbh
    
    ymax = (nbh * 4 * np.pi / 3)**(-1/3)

def Merger_Rate_Lite_wang(fbh, m, z):
    '''
    Get merger rate using wangsai's model, unit: 1/Gpc^3/yr
    '''
    
    Use_Cube = 0
    RhoCr = cosmo['RhoCr']
    zeq = cosmo['zeq']
    OmC = cosmo['OmC']
    H0 = cosmo['H0']
    msun = const['msun']
    c = const['c']
    G = const['G']
    
    RhoBH_eq = OmC * fbh * RhoCr * (1+zeq)**3
    mbh = m*msun # in kg
    nbh_eq = RhoBH_eq / mbh
    V = 1/nbh_eq # average volume for each BH
    if Use_Cube:
        x = V**(1/3)
    else:
        r = (3 * V/(4 * np.pi))**(1/3)
        x = 2*r
    
    T1 = (3/170) * c**5 * x**4
    T2 = (G * mbh)**3 * fbh**4
    T = T1/T2

    tc1 = (3/170) * c**5 * x**4 * fbh**(25/3)
    tc2 = (G * mbh)**3
    tc = tc1/tc2
    
    t = p21ct.z2t(z = z)
    
    if t < tc:
        dPdt = 3/58 * ((t/T)**(3/37) - (t/T)**(3/8))/t
    else:
        p1 = 3/58 * (t/T)**(3/8)
        p2 = fbh**(-29/8) * (t/tc)**(-29/56) - 1
        p3 = 1/t
        dPdt = p1 * p2 * p3
    
    R1 = 3 * H0**2 * fbh * OmC
    R2 = 8 * np.pi * G * mbh
    R3 = dPdt
    
    Gpc3yr = 9.274526196872258e+83
    
    R = R1*R2*R3 * Gpc3yr
    # R = R1*R2*R3

    return R
