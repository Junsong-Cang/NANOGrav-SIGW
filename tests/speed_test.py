from src.merger import *
import cProfile

fbh = 1e-2
mc = 1e-3
sbh = 0.3
v = np.logspace(0, 10, 100)
LineWidth = 2
FontSize = 18

from PyLab import *


def main():
    t1 = TimeNow()
    r3 = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh, 
        v = v,
        mf_model = 0, 
        sbh_width = 6, 
        nm = 50,
        nz = 50,
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
    return r3

if __name__ == '__main__':
    # Run the profiler
    cProfile.run('main()', sort='cumulative')
