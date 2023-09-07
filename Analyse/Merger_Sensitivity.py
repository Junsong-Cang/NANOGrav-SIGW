
LineWidth = 2
FontSize = 15

from PyLab import *
Curve_Path = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'

f1, h1 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.levitated_sensors.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f2, h2 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.BAW.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f3, h3 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.holometer.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f4, h4 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.EDGES.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f5, h5 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.ADMX.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f6, h6 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.SQMS.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f7, h7 = Read_Curve(
    File = Curve_Path + '2202.00695.fig1.ARCADE.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

Top = 1e-8
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.fill_between(f1, h1, Top,color = 'b', alpha=0.3, label='Levitated Sensors')
plt.fill_between(f2, h2, Top,color = 'r', alpha=0.3, label='BAW')
plt.fill_between(f3, h3, Top,color = 'g', alpha=0.3, label='holometer')
plt.fill_between(f4, h4, Top,color = 'y', alpha=0.3, label='EDGES')
plt.fill_between(f5, h5, Top,color = 'm', alpha=0.3, label='ADMX')
plt.fill_between(f6, h6, Top,color = 'c', alpha=0.3, label='SQMS')
plt.fill_between(f7, h7, Top,color = 'grey', alpha=0.3, label='ARCADE')

plt.xscale('log')
plt.yscale('log')

plt.xlabel('$f$ [MHz]',fontsize=FontSize,fontname='Times')
plt.ylabel('$h$ [strain]',fontsize=FontSize,fontname='Times')
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.title('Strain Sensitivities',fontsize=FontSize)
plt.ylim([1e-25, 1e-8])

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
