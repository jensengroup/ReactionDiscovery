%nprocshared=1
%mem=2GB
#opt=(calcall, ts, noeigen, nomicro, maxstep=15) 
external="/groups/kemi/mharris/github/xtb_ts_test/xtb_external.py"

something title

0 1
N -1.69996 1.18634 0.00248
C -0.57539 0.62929 -0.00039
N -0.42892 -0.80637 -0.00025
N -1.05913 -1.75125 -0.00301
O 1.37839 -0.86100 0.00341
N 1.70207 0.36445 -0.00003
C 0.73557 1.21867 -0.00293
H -2.51384 0.57601 0.00513
H 0.92715 2.27206 -0.00680

