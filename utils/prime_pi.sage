# number of special-q in [q0,q1] with factors in [qfac_min,qfac_max]
def nb_special_q(q0,q1,qfac_min=None,qfac_max=infinity,verbose=false):
   if qfac_min == None: # we consider only primes
      return prime_pi(q1)-prime_pi(q0)
   else:
      # determine maximum number of prime factors
      kmax = floor(log(q1)/log(qfac_min))
      assert kmax <= 3, "number of special-q factors is assumed <= 3"
      pi1 = prime_pi(min(q1,qfac_max))-prime_pi(max(q0,qfac_min))
      pi2 = 0
      pmax = min(qfac_max,floor(sqrt(q1)))
      for p in prime_range(qfac_min,pmax):
         pi2 += prime_pi(q1//p)-prime_pi(q0//p)
      pi3 = 0
      pmax = min(qfac_max,floor(q1^(1/3)))
      for p in prime_range(qfac_min,pmax):
         qmax = min(qfac_max,floor(sqrt(q1//p)))
         for q in prime_range(p,qmax):
            pi3 += prime_pi(q1//(p*q))-prime_pi(q0//(p*q))
      if verbose:
         print pi1, pi2, pi3
      return pi1 + pi2 + pi3
      
