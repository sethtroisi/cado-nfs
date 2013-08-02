# Sage functions to find optimize cofactorization strategies

# First you have to generate a set of methods with given B1,B2 bounds
# and possibly initial point (see README). This can be done by hand,
# or using the generate_methods() function below, for example:
# generate_methods("all.st", 500, 16, 18)

# Then, assuming this set of methods is in the file "all.st",
# the init() function will call the run_stat.sh program for each method,
# and produce corresponding output files st0, st1, ...
# Note: if a file, say st17, already exists, it is assumed to correspond to
# method of line 17 (starting at line 0) from all.st.

# Then the init_pr(bmin,bmax) function will set the initial probabilities.
# The arguments bmin and bmax specify the minimal and maximal bit size of
# the prime factors which are searched for, for example to search all prime
# factors of 29 to 33 bits:
# init_pr(29,33)

# Finally the find_best_chain(S, n, maxit) function will determine the optimal
# chain, with initial set of methods S, for a cafactor of n words, and up to
# maxit iterations, for example to find the best chain for one word with up to
# 30 methods (the global variable S is the initial set of methods):
# find_best_chain(S, 1, 30)

# read the file with given strategies
def init():
   global T, P, S, Method
   f = open ("all.st", "r")
   l = f.readlines()
   f.close()
   # launch run_stat.sh on each strategy
   Method = ["" for i in range(len(l))]
   for i in range(len(l)):
      Method[i] = l[i].strip('\n')
      try:
         f = open ("st" + str(i), "r")
      except:
         print "Generating data for method " + Method[i]
         command = "./run_stat.sh " + Method[i] + " > st" + str(i)
         os.system(command)
   # extract the timings and probabilities
   # T[i,n] is the time of strategy i for n words (1 <= n <= 3)
   T = dict()
   # P[i,t,k] is the probability of success of strategy i for t bits
   # (15 <= t <= 50) and k mod 12 (k = 1,5,7,11)
   P = dict()
   for i in range(len(l)):
      st = "st" + str(i)
      print "Parsing file " + st
      f = open(st, "r")
      li = f.readlines()
      for n in [1..3]:
	 found = 0
	 for s in li:
	    k = s.find(str(n) + " words")
	    if k <> -1:
	       found += 1
	       T[i,n] = float(s[k+9:-1])
	 if found <> 1:
	    print "Error, did not find exactly one occurrence of "+str(n)+" words"
	    raise ValueError
      for t in [15..50]:
	 found = 0
	 for s in li:
	    k = s.find(str(t) + ": ")
	    if k <> -1:
	       found += 1
	       P[i,t,1] = float(s[8:12])
	       P[i,t,5] = float(s[13:17])
	       P[i,t,7] = float(s[18:22])
	       P[i,t,11] = float(s[23:27])
	 if found <> 1:
	    print "Error, did not find exactly one occurrence of "+str(t)+":"
	    raise ValueError
      f.close()
   # initial set of methods
   S = Set(range(len(l)))

pr0 = dict()

# set initial probabilities
def init_pr(bmin,bmax):
   global pr, pr0
   assert 15 <= bmin <= bmax <= 50
   pr = dict()
   for t in [15..50]:
      for k in [1,5,7,11]:
	 pr0[t,k] = 1.0 / t / 4
         if t < bmin or bmax < t:
            pr0[t,k] = 0
         pr[t,k] = pr0[t,k]

# for p+1 with 2/7 as a starting point, we do
#     p-1 for p = 1, 7  mod 12
#     p+1 for p = 5, 11 mod 12

# for p+1 with 6/5 as a starting point, we do
#     p-1 for p = 1, 5  mod 12
#     p+1 for p = 7, 11 mod 12

def get_P(i,t,k):
   global pm1_last, pp1_27_last, pp1_65_last
   Pitk = P[i,t,k]
   Pitk_old = 0
   if Method[i][:3] == 'pm1' or (Method[i][:6] == 'pp1_27' and k in [1,7]) or (Method[i][:6] == 'pp1_65' and k in [1,5]):
      # check previous pm1, or pp1_27 if k=1,7, or pp1_65 if k=1,5
      if pm1_last >= 0:
         Pitk_old = P[pm1_last,t,k]
      if pp1_27_last >= 0 and k in [1,7]:
         Pitk_old = max(Pitk_old, P[pp1_27_last,t,k])
      if pp1_65_last >= 0 and k in [1,5]:
         Pitk_old = max(Pitk_old, P[pp1_65_last,t,k])
   if Method[i][:6] == 'pp1_27' and k in [5,11]:
      if pp1_27_last >= 0:
         Pitk_old = P[pp1_27_last,t,k]
      if pp1_65_last >= 0 and k in [7,11]:
         Pitk_old = max(Pitk_old, P[pp1_65_last,t,k])
   if Method[i][:6] == 'pp1_65' and k in [7,11]:
      if pp1_65_last >= 0:
         Pitk_old = P[pp1_65_last,t,k]
      if pp1_27_last >= 0 and k in [5,11]:
         Pitk_old = max(Pitk_old, P[pp1_27_last,t,k])
   return Pitk - Pitk_old

# find best method for n words
def find_best(S, n, verbose=True):
   global P, T, pr
   best_i = -1 # index of best method
   best_w = 0  # corresponding weight
   for i in S:
      # prob = 0
      failure = 1
      success = 0
      for t in [15..50]:
         for k in [1,5,7,11]:
            Pitk = get_P(i,t,k)
            if Pitk > 0:
               failure *= 1 - pr[t,k] * Pitk
      prob = 1 - failure
      w = -log(failure)/T[i,n]
      if verbose:
         print "method ", i, "proba:", prob, "cost:", T[i,n], "score:", w
      if w > best_w:
         best_i = i
         best_w = w
   return best_i, best_w

def update_pr(i):
   global P, pr
   for t in [15..50]:
      for k in [1,5,7,11]:
         # if P is the probability that it remains a prime of t bits
         # congruent to k mod 12, and p is the probability of success of
         # the method, then P becomes P*(1-p)/(1-P*p) after the method
         Pitk = get_P(i,t,k)
         pr[t,k] = pr[t,k]*(1-Pitk)/(1-pr[t,k]*Pitk)

def print_score():
   global pr, pr0
   print "remains: %.4f" % sum([pr[t,k] for t in [15..50] for k in [1,5,7,11]])
   for t in [15..50]:
      s0 = sum([pr0[t,k] for k in [1,5,7,11]])
      if s0 > 0:
         print "%1d:%.3f" % (t,sum([pr[t,k] for k in [1,5,7,11]])/s0),
   print

def find_best_chain(S, n, maxit, verbose=False):
   global Method, pm1_last, pp1_27_last, pp1_65_last
   pm1_last = -1
   pp1_27_last = -1
   pp1_65_last = -1
   print_score()
   total_time = 0
   for it in range(maxit):
      best_i, best_w = find_best (S, n, verbose)
      if best_i == -1:
         print "remaining methods do not give positive probability of success"
         for i in S:
            print Method[i]
         break
      print str(it) + ": best method " + str(best_i) + ": " + Method[best_i], "score:", float(best_w),
      update_pr(best_i)
      print_score()
      total_time += T[best_i,n]
      print "time so far:", total_time
      if Method[best_i][:3] == 'pm1':
         pm1_last = best_i
         S = S - Set([best_i])
      if Method[best_i][:6] == 'pp1_27':
         pp1_27_last = best_i
         S = S - Set([best_i])
      if Method[best_i][:6] == 'pp1_65':
         pp1_65_last = best_i
         S = S - Set([best_i])
      if Method[best_i][:3] == 'ecm' and Method[best_i][-3:] == ' 11':
         # we cannot replay the sigma=11 curves
         r = []
         for i in S:
            if Method[i][:3] == 'ecm' and Method[i][-3:] == ' 11':
               r.append(i)
         S = S - Set(r)

def estimate_chain(C, S, n, verbose=False):
   global Method, pm1_last, pp1_27_last, pp1_65_last
   pm1_last = -1
   pp1_27_last = -1
   pp1_65_last = -1
   print_score()
   total_time = 0
   for it in range(len(C)):
      best_i = C[it]
      if not (best_i in S):
         print "given method is not in remaining set:", best_i
         raise ValueError
      print str(it) + ": use given method", Method[best_i],
      update_pr(best_i)
      print_score()
      total_time += T[best_i,n]
      print "time so far:", total_time
      if Method[best_i][:3] == 'pm1':
         pm1_last = best_i
         S = S - Set([best_i])
      if Method[best_i][:6] == 'pp1_27':
         pp1_27_last = best_i
         S = S - Set([best_i])
      if Method[best_i][:6] == 'pp1_65':
         pp1_65_last = best_i
         S = S - Set([best_i])
      if Method[best_i][:3] == 'ecm' and Method[best_i][-3:] == ' 11':
         # we cannot replay the sigma=11 curves
         r = []
         for i in S:
            if Method[i][:3] == 'ecm' and Method[i][-3:] == ' 11':
               r.append(i)
         S = S - Set(r)

# Example: generate_methods("all.st", 500, 16, 18)         
def generate_methods(out, maxB1, mink, maxk):
   f = open (out, "w")
   B1 = 105
   while B1 < maxB1:
      oldB2 = 0
      for k in [mink..maxk]:
         B2 = k*B1
         k = ZZ(round(B2/210.))
         B2 = (2*k+1)*105
         if B2 > oldB2:
            oldB2 = B2
            f.write("pm1 " + str(B1) + " " + str(B2) + " 1\n")
            f.write("pp1_27 " + str(B1) + " " + str(B2) + " 1\n")
            f.write("pp1_65 " + str(B1) + " " + str(B2) + " 1\n")
            f.write("ecm " + str(B1) + " " + str(B2) + " 11\n")
            f.write("ecmm12 " + str(B1) + " " + str(B2) + " 2\n")
      B1 += B1//10
   f.close()
      


