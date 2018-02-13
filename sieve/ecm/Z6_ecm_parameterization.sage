
# Curve from Strafish on strike (theorem 5.4)
E1 = EllipticCurve([0,-1/2304,0,-5/221184,1/28311552])

# A curve is short Weierstrass form, isomorphic to E1 with integer coefficients
E2 = EllipticCurve([-9747, 285714])

# QQpqr.<p,q,r> = QQ[]
# I = QQpqr.ideal(E2.defining_polynomial()(z=r,x=p,y=q))

QQxyz.<x,y,z> = QQ[]
I = QQxyz.ideal(E2.defining_polynomial())

phi1 = E1.isomorphism_to(E2)
phi2 = E2.isomorphism_to(E1)

# Apply the Weierstrass transformation (u,r,s,t) in projective coordinate
# (x,y,z)  |-->  (U,V,W) = (u^2x + r*z , u^3y + su^2x + t*z, z).

U = (phi1).u^2*x+(phi1).r*z
V = (phi1).u^3*y+(phi1).s*(phi1).u^2*x+(phi1).t*z
W = z

l = lcm([c.denominator() for c in U.coefficients()\
         + V.coefficients() + W.coefficients()])

U *= l     # U = 144*(x + 3*z)
V *= l     # V = y
W *= l     # W = 2985984*z

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

eq_Ed = a*xE0^2 + yE0^2 - zE0^2 - d*tE0^2
print "P0 is on curve:       ",
print eq_Ed.numerator() in I\
    and not(eq_Ed.denominator() in I)\
    and (xE0*yE0 - zE0*tE0 in I)


# Compute Montgomery curve parameters
# following 20 years of ECM by Zimmermann and Dodson

xM0 = alpha3_n*beta3_d
zM0 = beta3_n*alpha3_d

#A = (((beta-alpha)^3)*(3*alpha+beta) / (4*alpha3*beta)) - 2
A = ((beta_n*alpha_d-alpha_n*beta_d)^3 * (3*alpha_n*beta_d + beta_n*alpha_d))\
     / (4*alpha3_n*beta_n*alpha_d*beta3_d) - 2
      
# Applying the morphism from the Montgomery curve By2 = x^3 + Ax^2 + x
# to its equivalent Edwards form ax^2 + y^2 = 1 + dx^2y^2
# sets a = (A+2)/B.
# We define B so that the corresponding Edwards curve has a = -1.
B = - (A+2)

# Check if PM0 = (xM0:yM0:zM0) is on the Montgomery curve equivalent to
# the Edwards curve eq_Ed
alpha3 = alpha^3
beta3 = beta^3
sigma = sigma_n / sigma_d
_b = alpha/beta3
y0 = (sigma^2-1)*(sigma^2-25)*(sigma^4-25)*alpha3_d*beta3_d

# Ad-hoc code for checking that b/B is a square in I
e2 = E2.defining_polynomial()
e3 = (-e2 + y^2*z)*z
f = (_b/(B*e3)).factor()
print "b/B is a square in I: ",\
    not(any([m[1] % 2 for m in list(f)])) and QQ(f.unit()).is_square()
yM02 = (_b/B)*y0^2
eq_M0 = B*yM02*zM0-xM0*(xM0^2+A*xM0*zM0+zM0^2)
print "M0 in on curve:       ",
print eq_M0.numerator() in I and not(eq_M0.denominator() in I)

# -----------------------------------------------------------------------------
# Define P3 = (xE3, yE3) a point of order 3 on eq_Ed
# Check that P3 in on the curve

xE3 = -r/((sigma-1)*(sigma+5))
yE3 = (alpha - beta)/(alpha + beta)
eq_P3 = a*xE3^2+yE3^2 - (1 + d*xE3^2*yE3^2)
print "P3 is on curve:       ",
print eq_P3.numerator() in I\
    and not(eq_P3.denominator() in I)

# Compute the coord of 3[P3] using tripling formula
def tpl () :
    X1 = xE3
    Y1 = yE3
    Z1 = 1
    YY = Y1^2
    aXX = a*X1^2
    Ap = YY+aXX
    B = 2*(2*Z1^2-Ap)
    xB = aXX*B
    yB = YY*B
    AA = Ap*(YY-aXX)
    F = AA-yB
    G = AA+xB
    xE = X1*(yB+AA)
    yH = Y1*(xB-AA)
    zF = Z1*F
    zG = Z1*G
    X3 = xE*zF
    Y3 = yH*zG
    Z3 = zF*zG
    T3 = xE*yH
    x3 = X3/Z3
    y3 = Y3/Z3
    eq_3P3 = a*x3^2+y3^2 - (1 + d*x3^2*y3^2)
    print "P3 is of order 3:     ",
    print eq_3P3.numerator() in I and not(eq_3P3.denominator() in I) 

tpl()
