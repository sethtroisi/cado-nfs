import sys

##### P+1
def pp1_parameterization (P0, p):
  F = GF(p)
  Fx = F['x']
  P = Fx([1, -P0, 1]) # P = x^2-P0*x+1
  if P.is_irreducible():
    F2 = F.extension (P, 'a')
    return F2, F2.gen()
  else:
    return F, P.roots()[0][0]


def Weierstrass_from_Montgomery (*args) :
    """
    Weierstrass_from_Montgomery (A, B)

    Weierstrass_from_Montgomery (A, B, [x, y])
    """
    if len (args) >= 2 :
        A = args[0]
        B = args[1]
        Ew = EllipticCurve ([0,A/B,0,1/B^2,0])
    if len (args) == 2:
        return Ew
    elif len (args) == 3 :
        x = args[2][0]
        y = args[2][1]
        return (Ew, Ew([x/B, y/B]))
    else :
        raise TypeError (sys._getframe().f_code.co_name + " takes 2 or 3 arguments")
    
# def MontgomeryCurve (*args):
#   """
#   MontgomeryCurve ([A,B])
#   MontgomeryCurve (Ring, [A,B])
#   """
#   if len(args) == 1:
#     A = args[0][0]
#     B = args[0][1]
#     return EllipticCurve ([0,A/B,0,1/B^2,0])
#   elif len(args) == 2:
#     R = args[0]
#     A = args[1][0]
#     B = args[1][1]
#     return EllipticCurve (R, [0,A/B,0,1/B^2,0])
#   else:
#     raise TypeError ("MontgomeryCurve takes at most 2 arguments")


def _Suyama_parameterization_internal (s):
  u = s^2-5
  v = 4*s
  x0 = u^3
  y0 = (s^2-1)*(s^2-25)*(s^4-25)
  z0 = v^3
  A = (v-u)^3*(3*u+v)/(4*u^3*v) - 2
  B = u/z0
  b = (v-u)^3*(3*u+v)/(16*u^3*v)                          # [ b = (A+2)/4 ]
  x3 = u
  y3 = 2*(s^2-25)*(s^2-1)*s/u
  z3 = v
  return A,B,b,x0,y0,z0,x3,y3,z3

def _Suyama_parameterization_check ():
  s = var('s')
  A,B,b,x0,y0,z0,x3,y3,z3 = _Suyama_parameterization_internal (s)
  assert (b-(A+2)/4).is_zero(), "b is not correct"
  assert (B*y0^2*z0 - (x0^3+A*x0^2*z0+x0*z0^2)).is_zero(), "P0 not on the curve"
  assert (B*y3^2*z3 - (x3^3+A*x3^2*z3+x3*z3^2)).is_zero(), "P3 not on the curve"
  assert (4*A*(x3/z3)^3+3*(x3/z3)^4+6*(x3/z3)^2-1).is_zero(), "not of order 3"

def Suyama_parameterization (*args, **kwargs):
  if len(args) == 1:
    s = args[0]
    curve_args = []
  elif len(args) == 2:
    s = args[1]
    curve_args = [args[0]]
  else:
    raise TypeError ("Suyama_parameterization takes at most 2 arguments")

  A,B,_,x0,y0,z0,x3,y3,z3 = _Suyama_parameterization_internal (s)
  curve_args.append([A,B])
  E = MontgomeryCurve (*curve_args)
  P0 = E(x0,y0,B*z0)
  if kwargs.get("check", False):
    P3 = E(x3,y3,B*z3)
    assert P3.order() == 3, "P3 is not of order 3"
  return E, P0

def _Z6_parameterization_internal (k):
  Eg = EllipticCurve (QQ, [-9747, 285714])
  Pg = Eg(15, 378)
  Qg = k*Pg
  p = Qg[0]/Qg[2]
  q = Qg[1]/Qg[2]

  sigma = -(p-213)/(p+3)

  alpha = (p-105)*(p+57)
  beta = q*(p-213)
  gamma = p^4 - 96*p^3 - 918*p^2 - 2808*p + 45346797
  delta = gamma + 2^5*3^7 *(p-33)^2
  epsilon =  p^2 - 66*p + 7569
  zeta = 2^2*3^6*q*(p-33)

  # Coeffs for "a=-1" Twisted Edwards curve
  a = -1
  d = -alpha^3*(p-51)*(p+327)/zeta^2
  xE0 = 2*3*(p+3)*alpha*beta*delta
  yE0 = (2*3^5)*(p-33)*alpha*gamma*epsilon
  zE0 = alpha^2*delta*epsilon
  tE0 = 2^2*3^6*(p-33)*(p+3)*beta*gamma
  xE3 = 2*3^2*q*alpha
  yE3 = 2*3^4*(p-33)*alpha
  zE3 = alpha^2
  tE3 = zeta
  
  # Coeffs for Montgomery curve
  A = 2*(a+d)/(a-d)
  B = 4/(a-d)
  b = zeta^2/(zeta^2-alpha^3*(p-51)*(p+327))               # [ b = (A+2)/4 ]
  xM0 = (p^2 + 114*p - 11331)^3
  yM0 = (alpha*epsilon*(epsilon+2^2*3^2*5*(p-105))^3)/(2*3*(p+3)*beta)
  zM0 = ((p-213)*(p+3))^3
  xM3 = 2*3^2*q*(epsilon+2^2*3^2*5*(p-105))
  yM3 = alpha*(epsilon+2^2*3^2*5*(p-105))
  zM3 = 2*3^2*(p+3)*beta

  Ed = (a, d, xE0, yE0, zE0, tE0, xE3, yE3, zE3, tE3)
  Mt = (A, B, b, xM0, yM0, zM0, xM3, yM3, zM3)
  return sigma, Ed, Mt



# Brent-Sumaya parameterization.
# Given rational parameter sigma, compute an elliptic curve in Weierstrass
# form together with a non torsion point.
def BrentSuyama_parameterization (sigma, p) :
    K = GF(p)
    u = K(sigma^2 - 5)
    v = K(4*sigma)
    a = (v - u)^3 * (3*u + v)
    b = 4*u^3*v
    A = a/b - 2
    # |A| != 2
    x = u^3 / v^3
    w = x^3 + A*x^2 + x
    # w != 0
    E = EllipticCurve (K, [0, w*A, 0, w^2, 0])
    return (E, E.point([x * w, w^2]))



# Montgomery elliptic parameterization with torsion group of order 12.
# Given integer parameter n, compute  an elliptic curve in Montgomery form
# with torsion subgroup of order 12.
# Returns an equivalent curve in short Weierstrass form.
def Monty12_parameterization (n, p, verbose=false) :
    K = GF(p)
    # v^2 = u^3 - 12u with generator (-2, 4)
    E = EllipticCurve (K, [0,0,0,-12,0])
    P0 = E.point([-2,4])
    [u,v,w] = n*P0
    if (verbose) :
        print u, ":", v
    t = v/(2*u)
    t2 = t^2
    a = (t2 - 1)/(t2 + 3)
    # a != 0
    X0 = 3*a^2+1
    Z0 = 4*a
    A = -(3*a^4+6*a^2-1)/(4*a^3)
    # |A| != 2
    B = (a^2-1)^2/(4*a^3)
    # B != 0

    # Montgomery curve is: B y^2 = x^3 + A x^2 + x
    # with y->B*Y, x->B*X we get Y^2 = X^3 + A/B*X^2 + 1/B^2*X
    x = X0/Z0
    y = 2 * (u^2+12)/(4*u)  / ((t^2+3)*(4*a));

    if (verbose) :
        print "t^2 = ", t2
        print "a =   ", a
        print "X0 =", X0, "is square (mod",p,")? ", X0.is_square()
        print "A = ", A, ", B = ", B
        print "x = ", x, ", y = ", y
        print Ew
        print "12 | order? ", Ew.order()%12 == 0
        # (a,a) is a point of order 3
        P3 = Ew([a*B,a*B^2])
        print "P3 = ", P3
        print "[3]P3 = ", 3*P3

    return Weierstrass_from_Montgomery (A, B, [x, y])



def Monty16_parameterization (n, p) :
    if n != 1:
        # Only one curve with torsion 16 implemented so far
        raise ValueError (sys._getframe().f_code.co_name + " requires n=1")
    K = GF(p)
    A = K(54721/14400)
    x = K(8/15)
    y = K(1)
    B = x^3 + A*x^2 + x
    return Weierstrass_from_Montgomery (A, B, [x, y])
        
        

def Twed12_parameterization (n,p, verbose=false) :
    K = GF(p)
    # Elliptic parameterization (parameter is called k).
    # E: y^2 = x^3 - 9747*x + 285714
    # rank = 1
    # P = (15, 378) is a point of infinite order
    # Torsion points: P2 = (33:0:1) of order 2
    #                 Q2 = (78:0:1) of order 2
    #                 P2+Q2 = (-111:0:1) of order 2
    # Reference: based on Theorem 5.4 of the article "Starfish on strike".
    # Valid parameter: n in Z \ { 0 }
    E = EllipticCurve (K, [-9747, 285714])
    P0 = E.point([15,378])
    x,y,z = n*P0
    
    U = 144*(x + 3*z)
    V = y
    W = 2985984*z
    
    sigma_d = 96*U                        # u0
    sigma_n = (W - sigma_d)               # u1
    sigma2_n = sigma_n^2                  # u2
    sigma2_d = sigma_d^2                  # u3
    alpha_n = sigma2_n - 5*sigma2_d       # u4
    alpha_d = sigma2_d                    # u3
    alpha3_n = alpha_n^3                  # u5
    alpha3_d = alpha_d^3                  # u6 = u3^3
    beta_n = 4*sigma_n                    # u7
    beta_d = sigma_d                      # u8 = u0
    beta3_n = beta_n^3                    # u9 = u7^3
    beta3_d = beta_d^3                    # u10 = u8^3
    gamma = (sigma_n - sigma_d)*(sigma_n + 5*sigma_d)*(sigma2_n + 5*sigma2_d)
    delta = (sigma2_n - 5*sigma2_d)^3
    epsilon = (beta_n * sigma_d)^3
    t0 = gamma*U^2
    t2 = 2*V*W*sigma_n*sigma_d^3

    zE0 = delta+epsilon
    xE0 = zE0
    zE0 = zE0*t0
    xE0 = xE0*t2
    yE0 = delta - epsilon
    tE0 = yE0
    yE0 = yE0*t0
    tE0 = tE0*t2

    # Check if PE0 = (xE0:yE0:zE0:tE0) is on the extended twisted
    # Edwards curve: (-1)*xE0^2 + yE0^2 - zE0^2 - d*tE0^2
    a = -1
    # Compute r = v/u^2 following Starfish on strike (Theorem 5.4)
    # with v = V/W and u = U/W
    r = (V*W)/U^2
    # Compute d following Starfish on strike (Theorem 5.4)
    alpha = alpha_n / alpha_d
    beta = beta_n / beta_d
    d = ((beta+alpha)^3*(beta-3*alpha)) / (r*(beta-alpha))^2

    # print "P0 is on curve:       ",
    # print (a*xE0^2 + yE0^2 - zE0^2 - d*tE0^2 == 0) and (xE0*yE0 == zE0*tE0)

    # Compute Montgomery curve parameters
    # following 20 years of ECM by Zimmermann and Dodson
    
    xM0 = alpha3_n*beta3_d
    zM0 = beta3_n*alpha3_d
    A = ((beta_n*alpha_d-alpha_n*beta_d)^3 * (3*alpha_n*beta_d + beta_n*alpha_d))\
        / (4*alpha3_n*beta_n*alpha_d*beta3_d) - 2
    
    # Applying the morphism from the Montgomery curve By2 = x^3 + Ax^2 + x
    # to its equivalent Edwards form ax^2 + y^2 = 1 + dx^2y^2
    # sets a = (A+2)/B.
    # We define B so that the corresponding Edwards curve has a = -1.
    B = -(A+2)
    
    # Check if PM0 = (xM0:yM0:zM0) is on the Montgomery curve equivalent to
    # the Edwards curve eq_Ed
    alpha3 = alpha^3
    beta3 = beta^3
    sigma = sigma_n / sigma_d
    _b = alpha/beta3
    y0 = (sigma^2-1)*(sigma^2-25)*(sigma^4-25)*alpha3_d*beta3_d
    
    # Ad-hoc code for checking that b/B is a square in I
#    e2 = E2.defining_polynomial()
#    e3 = (-e2 + y^2*z)*z
#    f = (_b/(B*e3)).factor()
#    print "b/B is a square in I: ",\
#        not(any([m[1] % 2 for m in list(f)])) and QQ(f.unit()).is_square()
    yM02 = (_b/B)*y0^2
#    eq_M0 = B*yM02*zM0-xM0*(xM0^2+A*xM0*zM0+zM0^2)
#    print "M0 in on curve:       ",
#    print eq_M0 == 0

    return Weierstrass_from_Montgomery (A, B, [xM0/zM0, K(sqrt(yM02)/zM0)])


def get_order_from_method(method, sigma, p) :
    # P-1
    if (method == 0) :
      R = GF(p)
      a = R(sigma)
      o = p-1
      po = a.multiplicative_order()
      lpf = largest_prime_factor(po)
      return (p, o, po, lpf)
    # P+1
    elif (method == 1):
      G, x0 = pp1_parameterization (sigma, p)
      o = p-1 if G == GF(p) else p+1
      po = x0.multiplicative_order()
      return (p, o, po, largest_prime_factor (po))
    # ECM - Brent-Suyama
    elif (method == 2) :
        T = BrentSuyama_parameterization (sigma, p)
    # ECM - Monty12
    elif (method == 3) :
        T = Monty12_parameterization (sigma, p)
    # ECM - Monty16
    elif (method == 4) :
        T = Monty16_parameterization (sigma, p)
    # ECM - Twed12
    elif (method == 5) :
        T = Twed12_parameterization (sigma, p)
        
    o = T[0].order()
    po = T[-1].order()
    lpf = largest_prime_factor (po)
    return (p, o, po, lpf)

    
def largest_prime_factor (n) :
    F = n.factor()
    return F[-1][0]
    
def is_powersmooth (fact, B):
  return all (p^e <= B for p,e in fact)

# Find the smallest prime p >= minp such that largest prime factor lpf
# in the order of the starting element for the factoring method satisfies
# minlpf <= lpf <= maxlpf, or if maxlpf == 0, minlpf <= lpf.
def prime_with_lpf_in_range (minp, minlpf, maxlpf, method, param):
  assert maxlpf == 0 or minlpf <= maxlpf, "maxlpf must be 0 or >= minlpf"
  p = next_prime (minp-1)
  while True:
    try:
      info = get_order_from_method (method, param, p)
      if minlpf <= info[3] and (maxlpf == 0 or info[3] <= maxlpf):
        return info
    except (ZeroDivisionError, ArithmeticError):
      pass
    p = next_prime (p)


# Find the smallest prime p >= minp such that the order of the starting element
# for the factoring method is B1,B2-powersmooth
def prime_with_B1_B2_smooth_order (minp, B1, B2, method, param):
  assert B2 > B1, "B2 must be > B1"
  p = next_prime (minp-1)
  while True:
    try:
      info = get_order_from_method (method, param, p)
      o = info[2]
      lpf = info[3]
      if B1 < lpf and lpf <= B2 and is_powersmooth(ZZ(o/lpf).factor(), B1):
        return info
    except (ZeroDivisionError, ArithmeticError):
      pass
    p = next_prime (p)
