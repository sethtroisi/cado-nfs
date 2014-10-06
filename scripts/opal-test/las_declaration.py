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
alambda = Parameter(kind='real', default=1.1, bound=(0.8, 1.4), name='alambda', description='Sieve survivor threshold')
rlambda = Parameter(kind='real', default=1.1, bound=(0.8, 1.4), name='rlambda', description='Sieve survivor threshold')
alim = Parameter(kind='integer', default=600000, bound=(50000, 2000000), name='alim', description='Factor base bound, algebraic side')
rlim = Parameter(kind='integer', default=600000, bound=(30000, 2000000), name='rlim', description='Factor base bound, rational side')
lpba = Parameter(kind='integer', default=22, bound=(18, 25), name='lpba', description='Large prime bound, algebraic side')
lpbr = Parameter(kind='integer', default=22, bound=(18, 25), name='lpbr', description='Large prime bound, rational side')
mfba = Parameter(kind='integer', default=25, bound=(20, 30), name='mfba', description='Cofactorization bound, algebraic side')
mfbr = Parameter(kind='integer', default=25, bound=(20, 30), name='mfbr', description='Cofactorization bound, rational side')
I = Parameter(kind='integer', default=11, bound=(10, 12), name='I', description='Sieve region size')
LAS.add_param(alambda)
LAS.add_param(rlambda)
LAS.add_param(alim)
LAS.add_param(rlim)
LAS.add_param(lpba)
LAS.add_param(lpbr)
LAS.add_param(mfba)
LAS.add_param(mfbr)
LAS.add_param(I)

# Define relevant measure and register with algorithm.
sievetime = Measure(kind='real', name='SIEVETIME', description='Time in the sieving')
rels = Measure(kind='integer', name='RELATIONS', description='Relations found in the sieving')
LAS.add_measure(sievetime)
LAS.add_measure(rels)
