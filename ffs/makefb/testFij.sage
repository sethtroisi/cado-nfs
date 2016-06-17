load "readparam.sage"

def skewgauss(q, rho, S):
    Fpt = q.parent()
    t = Fpt.gen()
    a0 = q
    b0 = Fpt(0)
    a1 = rho
    b1 = t^S
    while True:
        qq = a0 // a1
        if qq.degree()+b1.degree() > a0.degree():
            break
        a0 -= qq*a1
        b0 -= qq*b1
        qq = a1 // a0
        if qq.degree()+b0.degree() > a1.degree():
            break
        a1 -= qq*a0
        b1 -= qq*b0
        if a0.degree() <= b0.degree():
            break
    return a0, b0 // t^S, a1, b1 // t^S

# this will define, among others, pol0
readparam("../param.2.607")
t = q0.parent().gen()
x = pol0.parent().gen()

q = t^25 + t^3 + 1
rho = t^24 + t^23 + t^22 + t^21 + t^19 + t^16 + t^15 + t^13 + t^11 + t^10 + t^7 + t^3

a0, b0, a1, b1 = skewgauss(q, rho, 0)

a = x*a0+a1
b = x*b0+b1

Fij = pol0.parent()(0)
C = pol0.coeffs()
d = len(C)-1
for k in range(0, len(C)):
    Fij += C[k]*a^k*b^(d-k)

def int2pol(i):
   l = ZZ(i).bits()
   T = t.parent()
   return T(sum([l[k]*t^k for k in range(len(l))]))

# check for cancellations in the evaluation of Fij for i < 2^I, j < 2^J
# return the number of times each coefficients gives the degree of the sum
# example: check_cancellation(Fij,9,9)
def check_cancellation(Fij, I, J, verbose=False):
   fij = Fij.coefficients()
   one = t//t
   # precompute polynomials
   L = [int2pol(i) for i in range(1,2^max(I,J))]
   maxi = [0 for k in range(d+1)]
   for ti in L:
      # precompute fij[k]*ti^k
      a = [fij[k]*ti^k for k in range(d+1)]
      for tj in L:
	 if ti.gcd(tj) == one:
	    dij = [a[k]*tj^(d-k) for k in range(d+1)]
	    dijs = [dij[k].degree() for k in range(d+1)]
	    s = sum(dij)
            m = max(dijs)
            if verbose:
               print dijs, s.degree()
            for k in range(d+1):
               if dijs[k] == m:
                  maxi[k] += 1
	    if s.degree() <> m:
	       print "cancellation for i=", i, "j=", j, m - s.degree()
   return maxi
