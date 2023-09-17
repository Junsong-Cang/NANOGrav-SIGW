from src.merger import *

def Get_dEGW_dv_2(v, eta, M):
    
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
    
    print('v1 = ', v1, ', v2 = ', v2, ', v3 = ', v3, ', v4 = ', v4, ', M = ', M/2)

    r = (np.pi * G)**(2/3) * (M*msun)**(5/3) * eta * F
    return r

m = 1
v = 10
Get_dEGW_dv_2(v = 132, eta = 1/4, M=2)
v1 =  4050.381434044801

def test_ergf(z):
    
    vr = v * (1+z)
    dif = vr - v1
    return dif


z1 = PyLab.Solve(F = test_ergf, Xmin = 0, Xmax = 1000, Precision=1e-3, show_status=1)
print(z1)
