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
lim0 = Parameter(kind='integer', default=600000, bound=(30000, 2000000), name='lim0', description='Factor base bound, rational side')
lim1 = Parameter(kind='integer', default=600000, bound=(50000, 2000000), name='lim1', description='Factor base bound, algebraic side')
lpb0 = Parameter(kind='integer', default=22, bound=(18, 25), name='lpb0', description='Large prime bound, rational side')
lpb1 = Parameter(kind='integer', default=22, bound=(18, 25), name='lpb1', description='Large prime bound, algebraic side')
mfb0 = Parameter(kind='integer', default=25, bound=(20, 30), name='mfb0', description='Cofactorization bound, rational side')
mfb1 = Parameter(kind='integer', default=25, bound=(20, 30), name='mfb1', description='Cofactorization bound, algebraic side')
lambda0 = Parameter(kind='real', default=1.1, bound=(0.8, 1.4), name='lambda0', description='Sieve survivor threshold')
lambda1 = Parameter(kind='real', default=1.1, bound=(0.8, 1.4), name='lambda1', description='Sieve survivor threshold')
ncurves0 = Parameter(kind='integer', default=6, bound=(3,10), name='ncurves0', description='Cofactorization curves, side 0')
ncurves1 = Parameter(kind='integer', default=6, bound=(3,10), name='ncurves1', description='Cofactorization curves, side 1')
I = Parameter(kind='integer', default=11, bound=(10, 12), name='I', description='Sieve region size')
LAS.add_param(lim0)
LAS.add_param(lim1)
LAS.add_param(lpb0)
LAS.add_param(lpb1)
LAS.add_param(mfb0)
LAS.add_param(mfb1)
LAS.add_param(lambda0)
LAS.add_param(lambda1)
LAS.add_param(ncurves0)
LAS.add_param(ncurves1)
LAS.add_param(I)

# Define relevant measure and register with algorithm.
sievetime = Measure(kind='real', name='SIEVETIME', description='Time in the sieving')
rels = Measure(kind='integer', name='RELATIONS', description='Relations found in the sieving')
LAS.add_measure(sievetime)
LAS.add_measure(rels)
