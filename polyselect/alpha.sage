
plimit_rootprops=1000
# primes_rootprops=prime_range(2,200)

# Also counts projective roots. The input must be a polynomial.
def allroots(P,p):
    Pp= PolynomialRing(GF(p),'x')(P)
    s=sum([r[1] for r in Pp.roots()])
    s+=P.degree()-Pp.degree()
    return s

def val0(P0,p):
   v = valuation (P0.content(), p)
   ZP= P0.parent()
   x = ZP.gen()
   P = P0.parent()(P0/p^v)
   Q = GF(p)['x'](P.derivative())
   for r in P.roots(GF(p)):
       if Q(r[0]) != 0: v += 1/(p-1)
       else:
         P2 = P(Integers()(r[0])+p*x)
         v += val0(P2, p)/p
   return v

def val(P,p):
    r = val0(P, p) * p
    ZP= P.parent()
    x = ZP.gen()
    r += val0((P.reverse())(p*x), p)
    r /= p+1
    return r

def get_alpha(f,B):
   s = 0
   disc = f.discriminant()
   for p in prime_range(2,B):
     if disc % p == 0:
         # e := est_valuation (f, p, 1, 99, 1, 99);
         e = val(f,p)
         ds = float((1/(p-1)-e)*log(p))
#        lprint(p," divides disc(f)", evalf(e,3));
     else:
         q = 0
         q += allroots(f,p)
         ds = float((1-q*p/(p+1))*log(p)/(p-1))
     # print "alpha(%d)=%f" % (p,ds)
     s += ds

   return s

