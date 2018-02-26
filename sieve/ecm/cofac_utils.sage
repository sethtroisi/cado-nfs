
def Weierstrass_from_Montgomery (A, B, x, y) :
    Ew = EllipticCurve ([0,A/B,0,1/B^2,0])
    xw = x/B
    yw = y/B
    assert (Ew.is_on_curve(xw,yw))
    return [Ew, Ew([xw,yw])]




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
    assert (A != 2 and A != -2)
    x = u^3 / v^3
    w = x^3 + A*x^2 + x
    assert (w != 0)
    E = EllipticCurve (K, [0, w*A, 0, w^2, 0])
    return [E, E.point([x * w, w^2])];



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
    assert (t2 != -3)
    a = (t2 - 1)/(t2 + 3)
    assert (a != 0)
    X0 = 3*a^2+1
    Z0 = 4*a
    A = -(3*a^4+6*a^2-1)/(4*a^3)
    assert (A != 2 or A != -2)
    B = (a^2-1)^2/(4*a^3)
    assert (B != 0)

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

    # Ew = EllipticCurve(K, [0,A/B,0,1/B^2,0])
    # xw = x/B
    # yw = y/B
    # assert (Ew.is_on_curve(xw, yw))
    
        
    return Weierstrass_from_Montgomery (A, B, x, y)



# def Monty16_parameterization (n, p) :
# TODO



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

    print "P0 is on curve:       ",
    print (a*xE0^2 + yE0^2 - zE0^2 - d*tE0^2 == 0) and (xE0*yE0 == zE0*tE0)

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
    eq_M0 = B*yM02*zM0-xM0*(xM0^2+A*xM0*zM0+zM0^2)
    print "M0 in on curve:       ",
    print eq_M0 == 0

    return Weierstrass_from_Montgomery (A, B, xM0/zM0, K(sqrt(yM02)/zM0))

            

    
# get_order

# P-1
# Inputs: p, param
# Output: [p, o, po, lpf]

# R = IntegerModRing (p-1)
# a = R(param)
# o = p-1
# po = a.order()
# lpf = largest_prime_factor(po)
# return [p, o, po, lpf]

def get_order_from_method(method, sigma, p) :
# ECM
# Inputs:
# Output: [p, o, po, lpf]

# Brent-Suyama (sigma)
# Compute elliptic curve E parameterized by sigma and base point (x,y)
    if (method == 2) :
        T = BrentSuyama_parameterization (sigma, p)
    elif (method == 3) :
        T = Monty12_parameterization (sigma, p)
        #elif (method == 4) :
        #    T = Monty16_parameteriztion (n, p)
    elif (method == 5) :
        T = Twed12_parameterization (sigma, p)
        
        # T = [Elliptic curve E in Weierstrass form, curve parameter, non torsion point on E]
        
    o = T[0].order()
    po = T[-1].order()
    lpf = largest_prime_factor (po)
    return p, o, po, lpf

    
def largest_prime_factor (n) :
    F = n.factor()
    F.sort()
    return F[-1][0]
    
    
