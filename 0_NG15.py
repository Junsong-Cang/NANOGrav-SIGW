# Sample NG15 data with no additional priors

FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/0_NG15/0_NG15_'
p1_info = {'name':'LgF', 'min':-11, 'max':-5, 'latex':'\log_{10} (f_*/{\mathrm{Hz}})'}
p2_info = {'name':'LgA', 'min':-3, 'max':1, 'latex':'\log_{10} A'}
p3_info = {'name':'S', 'min':0.1, 'max':3, 'latex':'\Delta'}

ndim = 3

from src.main import *
from pymultinest.solve import solve
from PyLab import *
do_mcmc = eval(input('Do you want to do mcmc?\n'))

def Flat_Prior(cube):
  # A flat prior
	p1, p2, p3 = cube
	p1 = p1*(p1_info['max']-p1_info['min']) + p1_info['min']
	p2 = p2*(p2_info['max']-p2_info['min']) + p2_info['min']
	p3 = p3*(p3_info['max']-p3_info['min']) + p3_info['min']
	NewCube = [p1, p2, p3]
	return NewCube

def LogLike(theta):
  LgF, LgA, S = theta
  LnL = Log_Like(
    LgFc = LgF,
    LgA = LgA,
    Sigma = S,
    GW_Data_Set = [1, 0, 0, 0, 0, 0],
    Use_Optimistic = 0,
    Use_Logspace = True,
    Fig_Was_2_Sigma = False,
    Use_Neff = False,
    Use_Fbh = False,
    chi2_format = 0)
  return LnL

if do_mcmc:
  # run MultiNest
  t1 = TimeNow()
  r = solve(
      LogLikelihood=LogLike,
      Prior=Flat_Prior,
	    n_dims=ndim,
      resume = False,
      outputfiles_basename=FileRoot, 
      verbose=True,
      n_iter_before_update = 10)
  Timer(t1)
else:
  # Print mcmc info
  info = [p1_info, p2_info, p3_info]
  print_mcmc_info(FileRoot, info)
  Add_derived_param(FileRoot = FileRoot)
