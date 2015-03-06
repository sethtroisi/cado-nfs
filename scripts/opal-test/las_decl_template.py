# Description of the foward finite-difference "algorithm".
from opal.core.algorithm import Algorithm
from opal.core.parameter import Parameter
from opal.core.measure   import Measure

# Define Algorithm object.
LAS = Algorithm(name='LAS', description='Lattice Siever')

# Register executable for LAS.
LAS.set_executable_command('python las_run.py')

# Register parameter file used by black-box solver to communicate with LAS.
#LAS.set_parameter_file('las.param')  # Should be chosen automatically and hidden.

# Define parameter and register it with algorithm.
rlim = Parameter(kind='integer', default=rlim_def, bound=(rlim_min, rlim_max), name='rlim', description='Factor base bound, rational side')
alim = Parameter(kind='integer', default=alim_def, bound=(alim_min, alim_max), name='alim', description='Factor base bound, algebraic side')
lpbr = Parameter(kind='integer', default=lpbr_def, bound=(lpbr_min, lpbr_max), name='lpbr', description='Large prime bound, rational side')
lpba = Parameter(kind='integer', default=lpba_def, bound=(lpba_min, lpba_max), name='lpba', description='Large prime bound, algebraic side')
mfbr = Parameter(kind='integer', default=mfbr_def, bound=(mfbr_min, mfbr_max), name='mfbr', description='Cofactorization bound, rational side')
mfba = Parameter(kind='integer', default=mfba_def, bound=(mfba_min, mfba_max), name='mfba', description='Cofactorization bound, algebraic side')
ncurves0 = Parameter(kind='integer', default=ncurves0_def, bound=(ncurves0_min,ncurves0_max), name='ncurves0', description='Cofactorization curves, side 0')
ncurves1 = Parameter(kind='integer', default=ncurves1_def, bound=(ncurves1_min,ncurves1_max), name='ncurves1', description='Cofactorization curves, side 1')
I = Parameter(kind='integer', default=I_def, bound=(I_min, I_max), name='I', description='Sieve region size')
LAS.add_param(rlim)
LAS.add_param(alim)
LAS.add_param(lpbr)
LAS.add_param(lpba)
LAS.add_param(mfbr)
LAS.add_param(mfba)
LAS.add_param(ncurves0)
LAS.add_param(ncurves1)
LAS.add_param(I)

# Define relevant measure and register with algorithm.
sievetime = Measure(kind='real', name='SIEVETIME', description='Time in the sieving')
rels = Measure(kind='integer', name='RELATIONS', description='Relations found in the sieving')
LAS.add_measure(sievetime)
LAS.add_measure(rels)
