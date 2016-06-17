
def FindGroupOrder(p,s):
   K = GF(p)
   v = K(4*s)
   u = K(s**2-5)
   x = u**3
   b = 4*x*v
   a = (v-u)**3*(3*u+v)
   A = a/b-2
   x = x/v**3
   b = x**3 + A*x**2 + x
   E = EllipticCurve(K,[0,b*A,0,b**2,0])
   return factor(E.cardinality())

def is_smooth(n,B1,B2,exp2):
   l = factor(n)
   m = len(l)
   for i in range(m-1):
      if l[i][0] == 2:
         extra = 2^exp2
      else:
         extra = 1
      if l[i][0]^l[i][1] > extra*B1:
         return False
   if l[m-1][0]^l[m-1][1] > B2:
         return False
   return True

def pm1(p,B1,B2):
   return is_smooth(p-1,B1,B2,0)

def pp1(p,B1,B2):
   x0 = 2/7
   delta = (x0^2-4) % p
   return is_smooth(p-delta.jacobi(p),B1,B2,0)

def brent12(p,B1,B2,sigma):
   l = FindGroupOrder(p,sigma)
   return is_smooth(l.expand(),B1,B2,2)

def monty12(p,B1,B2,k,extra):
   E=EllipticCurve([-12,0])
   P=E(-2,4)
   kP=k*P; u=kP[0]; v=kP[1]
   t2 = (u^2-12)/(4*u)
   a = (u^2 - 4*u - 12)/(u^2 + 12*u - 12)
   A = (-3*a^4-6*a^2+1)/(4*a^3)
   B = (a^2-1)^2/(4*a^3)
   x = (3*a^2+1)/(4*a)
   # curve is B y^2 = x^3 + A x^2 + x
   # with y->B*Y, x->B*x we get Y^2 = X^3 + A/B*X^2 + 1/B^2*X
   E = EllipticCurve(GF(p), [0,A/B,0,1/B^2,0])
   l = E.cardinality()
   return is_smooth(l,B1,B2,extra)

# return True iff p is found by facul with NB_CURVES=n
def cofactor_prime(p, n):
   if pm1(p,315,2205):
      return True
   if pp1(p,525,3255):
      return True
   if monty12(p,105,3255,2,2):
      return True
   if n > 0 and brent12(p,315,5355,11):
      return True
   B1 = 105.0
   for i in range(4,n+3):
      B1 += sqrt (B1)
      B2 = 17.0 * B1
      k = ZZ((B2 / 210.0).trunc())
      if monty12(p,ZZ(B1.trunc()), (2*k+1)*105, i-1, 1):
         return True
   return False

def cofactor(q, n):
   l = factor(q)
   for p,e in l:
      if cofactor_prime(p, n):
         return True
   return False

# returns probability of finding a prime p < 2^lpb with NB_CURVES=n
# on K largest primes below 2^lpb
# lpb n:proba
# 22  0:0.7560 1:0.9074
# 23  1:0.8591 2:0.9059
# 24  2:0.8572 3:0.8990
# 25  4:0.8878 5:0.9194
# 26  5:0.8715 6:0.9065
# 27  7:0.8777 8:0.9053
# 28  9:0.8743 10:0.9010
# 29  12:0.8863 13:0.9091
# 30  15:0.8939 16:0.9134
# 31  17:0.8866 18:0.9039
# 32  20:0.8893 21:0.9076
# 33  24:0.8963 25:0.9121
def proba_cofactor(lpb, n, K=10000):
   s = 0
   p = ZZ(round(2^(lpb-0.5)))
   for k in range(K):
      p = next_prime (p)
      if cofactor_prime(p, n):
         s += 1
   return 1.0*s/K
   

def do_table(lpbmin, lpbmax, **kwargs):
    ncurves=0
    ntries=100
    target_prob=0.9
    if "ncurves" in kwargs:
        ncurves=kwargs["ncurves"]
    if "ntries" in kwargs:
        ntries=kwargs["ntries"]
    if "target_prob" in kwargs:
        target_prob=kwargs["target_prob"]
    for lpb in range(lpbmin, lpbmax+1):
        results=[]
        printed=[]
        ncurves = ncurves - 1
        while len(results) < 1 or results[-1][1] < target_prob:
            ncurves = ncurves + 1
            p=proba_cofactor(lpb, ncurves, ntries)
            # print("\t\t/* %d,%d:%f */" % (lpb, ncurves, p))
            printed.append("%d:%f" % (ncurves, p))
            results.append((ncurves, p))
        if len(printed) > 2:
            printed=printed[-2:]
        print("/* lpb=%d */ %d, /* %s */" % (lpb, ncurves, ", ".join(printed)))

do_table(10,64,ntries=10000)


