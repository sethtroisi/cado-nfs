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

    # The curve BY^2Z = X(X^2 + AXZ + Z^2) has a torsion subgroup of order 12

    x = X0/Z0
    y = 2 * (u^2+12)/(4*u)  / ((t^2+3)*(4*a));
    Ew = EllipticCurve (K, [0,A*B,0,B^2, 0])
    xw = B*x
    yw = B^2*y
    assert (Ew.is_on_curve(xw, yw))

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
        
    return [Ew, Ew([xw,yw])]

# def Monty16_parameterization (n, p) :

# def Twed12_parameterization (n,p) :
    

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
        #elif (method == 5) :
        #    T = Twed12_parameteriation (n, p)
        
        # T = [Elliptic curve E in Weierstrass form, curve parameter, non torsion point on E]
        
    o = T[0].order()
    po = T[-1].order()
    lpf = largest_prime_factor (po)
    return [p, o, po, lpf]

    
def largest_prime_factor (n) :
    F = n.factor()    F.sort()
    return F[-1][0]
    
    
