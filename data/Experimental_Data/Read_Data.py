LineWidth = 2
FontSize = 18

import numpy as np
import matplotlib.pyplot as plt

def Get_Data(anchor_file, x_axis_file, top_file, middle_file, lower_file):
    '''
    Read data and get error bars and mediam
    '''
    def pixel_2_real(xp1, xp2, xr1, xr2, xp):
        # go from pixel to real coordinate
        xr = (xr2 - xr1)*(xp - xp1)/(xp2 - xp1) + xr1
        return xr

    anchor = np.loadtxt(anchor_file)
    x_axis_pixel = np.loadtxt(x_axis_file)
    top_pixel = np.loadtxt(top_file)
    lower_pixel = np.loadtxt(lower_file)
    middle_pixel = np.loadtxt(middle_file)
    
    # Load x-axis
    x1r = anchor[0]
    x1p = anchor[2]
    x2r = anchor[4]
    x2p = anchor[6]
    x_axis = pixel_2_real(xp1 = x1p, xp2 = x2p, xr1 = x1r, xr2 = x2r, xp = x_axis_pixel)

    # load y
    y1r = anchor[9]
    y1p = anchor[11]
    y2r = anchor[13]
    y2p = anchor[15]
    
    # top
    top = pixel_2_real(xp1 = y1p, xp2 = y2p, xr1 = y1r, xr2 = y2r, xp = top_pixel)
    # mid
    mid = pixel_2_real(xp1 = y1p, xp2 = y2p, xr1 = y1r, xr2 = y2r, xp = middle_pixel)
    # lower
    lower = pixel_2_real(xp1 = y1p, xp2 = y2p, xr1 = y1r, xr2 = y2r, xp = lower_pixel)
    
    result = {'x_axis' : x_axis, 'middle': mid, 'top' : top, 'lower' : lower}
    
    return result

# Get NG15 data

NG15_conservative = Get_Data(
    anchor_file = '2306.16215_NG15_anchor.txt',
    x_axis_file = '2306.16215_NG15_x_axis.txt',
    top_file = '2306.16215_NG15_top_conservative.txt',
    middle_file = '2306.16215_NG15_middle.txt',
    lower_file = '2306.16215_NG15_lower_conservative.txt')

NG15_optimistic = Get_Data(
    anchor_file = '2306.16215_NG15_anchor.txt',
    x_axis_file = '2306.16215_NG15_x_axis.txt',
    top_file = '2306.16215_NG15_top_optimistic.txt',
    middle_file = '2306.16215_NG15_middle.txt',
    lower_file = '2306.16215_NG15_lower_optimistic.txt')

IPTA_optimistic = Get_Data(
    anchor_file = '2306.16215_IPTA_anchor.txt',
    x_axis_file = '2306.16215_IPTA_x_axis.txt',
    top_file = '2306.16215_IPTA_top_optimistic.txt',
    middle_file = '2306.16215_IPTA_middle.txt',
    lower_file = '2306.16215_IPTA_lower_optimistic.txt')

IPTA_conservative = Get_Data(
    anchor_file = '2306.16215_IPTA_anchor.txt',
    x_axis_file = '2306.16215_IPTA_x_axis.txt',
    top_file = '2306.16215_IPTA_top_conservative.txt',
    middle_file = '2306.16215_IPTA_middle.txt',
    lower_file = '2306.16215_IPTA_lower_conservative.txt')

PPTA = Get_Data(
    anchor_file = '2306.17124_anchor.txt',
    x_axis_file = '2306.17124_PPTA_x_axis.txt',
    top_file = '2306.17124_PPTA_top.txt',
    middle_file = '2306.17124_PPTA_middle.txt',
    lower_file = '2306.17124_PPTA_lower.txt')

EPTA = Get_Data(
    anchor_file = '2306.17124_anchor.txt',
    x_axis_file = '2306.17124_EPTA_x_axis.txt',
    top_file = '2306.17124_EPTA_top.txt',
    middle_file = '2306.17124_EPTA_middle.txt',
    lower_file = '2306.17124_EPTA_lower.txt')

data = NG15_conservative
np.savez('NG15_conservative.npz', 
         freq = 10**data['x_axis'], 
         LgGW_middle = data['middle'],
         LgGW_top = data['top'], 
         LgGW_lower = data['lower'])

data = NG15_optimistic
np.savez('NG15_optimistic.npz', 
         freq = 10**data['x_axis'], 
         LgGW_middle = data['middle'], 
         LgGW_top = data['top'], 
         LgGW_lower = data['lower'])

data = IPTA_optimistic
np.savez('IPTA_optimistic.npz', 
         freq = 10**data['x_axis'], 
         LgGW_middle = data['middle'], 
         LgGW_top = data['top'], 
         LgGW_lower = data['lower'])

data = IPTA_conservative
np.savez('IPTA_conservative.npz', 
         freq = 10**data['x_axis'], 
         LgGW_middle = data['middle'], 
         LgGW_top = data['top'], 
         LgGW_lower = data['lower'])

data = PPTA
np.savez('PPTA.npz', 
         freq = 10**data['x_axis'], 
         LgGW_middle = data['middle'], 
         LgGW_top = data['top'], 
         LgGW_lower = data['lower'])

data = EPTA
np.savez('EPTA.npz', 
         freq = 10**data['x_axis'], 
         LgGW_middle = data['middle'], 
         LgGW_top = data['top'], 
         LgGW_lower = data['lower'])

'''
# Get some plots
r1 = NG15_conservative
r2 = NG15_optimistic
r1 = IPTA_conservative
r2 = IPTA_optimistic

f = 10**r1['x_axis']
m = r1['middle']
t1 = r1['top']
l1 = r1['lower']
st1 = t1 - m
sl1 = m - l1

t2 = r2['top']
l2 = r2['lower']
st2 = t2 - m
sl2 = m - l2

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
plt.errorbar(f, m, [sl1, st1], color = 'k', linewidth=LineWidth, label = 'Conservative', fmt='+')
#plt.errorbar(f, m, [sl2, st2], color = 'r', linewidth=LineWidth, label = 'Optimistic', fmt='+')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$T$ [K]',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.tight_layout()

plt.show()
'''


# Get some plots
r = EPTA
f = 10**r['x_axis']
m = r['middle']
t = r['top']
l = r['lower']
st = t - m
sl = m - l

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
plt.errorbar(f, m, [sl, st], color = 'k', linewidth=LineWidth, label = 'Conservative', fmt='+')
#plt.errorbar(f, m, [sl2, st2], color = 'r', linewidth=LineWidth, label = 'Optimistic', fmt='+')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$T$ [K]',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.tight_layout()

plt.show()
