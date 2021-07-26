%nprocshared=1
%mem=2GB
#opt=(calcall, ts, noeigen, nomicro, maxstep=15) 
external="/groups/kemi/mharris/github/xtb_ts_test/xtb_external.py"

something title

0 1
N -0.513432 0.462212 -0.750430
C -0.343250 -0.576329 -0.099047
N -1.243208 -1.646316 0.037244
N -1.678034 -2.433563 0.698518
O 1.081206 -2.061815 -1.041629
N 1.588715 -1.661934 0.076827
C 0.894743 -0.851923 0.728812
H -1.291932 0.542192 -1.376505
H 1.056510 -0.504695 1.757240

