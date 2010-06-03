load "cado-nfs/trunk/polyselect/alpha.sage"

##############################################################
#TQponctuel                                                  #
##############################################################
#Input:                                                      #
#  - an integer n                                            #
#  - a prime integer p.                                      #
#We assume that n is a quadratic residue (mod p)             #
#                                                            #
#Output:                                                     #
#  - two pairs of quadratic polynomials defined from n and p #
#For each pair (f,g), f and g have a common root modulo n    #
##############################################################

def TQponctuel(n, p):
    
    z = ZZ(isqrt(n))
    listpoly = []
    
    #We compute the vectors a, b, c the way Peter Montgomery explained, starting from a root of n modulo p
    #These vectors span a lattice of Z^3

    c0 = p
    c1 = ZZ(mod(n,p).sqrt())
    c1 += p * ( ZZ(z + p//2 - c1)//p )
    c2 = ((c1^2) - n)//p
    c = [c0, c1, c2]
    
    a = [c1, -p, 0]
    
    s = ZZ( 1/ (mod(c1, p) ) )
    b0 = (c1*((c2*s)%p) - c2) // p
    b1 = - ((c2*s)%p)
    b = [b0, b1, 1]
    

    #We then perform the LLL reduction on a and b
    A = Matrix([[a[0], a[1], a[2]], [b[0], b[1], b[2]]])
    B = A.LLL()
    
    R.<x> = PolynomialRing(ZZ)
    f = B[0,0] + B[0,1]*x + B[0,2]*x^2
    g = B[1,0] + B[1,1]*x + B[1,2]*x^2
    
    #If the resultant of the polynomials is n, the polynoms are valid
    #I wonder if this last verification is useful, it seems bo always be the case
    res = f.resultant(g)
    if res == n:
        listpoly.append([f,g])


    #We start all over again, taking this time the other square root of n modulo p
    c1 = -c1
    c1 += p * ( ZZ(z + p//2 - c1)//p )
    c2 = ((c1^2) - n)//p
    c = [c0, c1, c2]
    
    a = [c1, -p, 0]
    
    s = ZZ( 1/ (mod(c1, p) ) )
    b0 = (c1*((c2*s)%p) - c2) // p
    b1 = - ((c2*s)%p)
    b = [b0, b1, 1]
    
    A = Matrix([[a[0], a[1], a[2]], [b[0], b[1], b[2]]])
    B = A.LLL()
    
    f = B[0,0] + B[0,1]*x + B[0,2]*x^2
    g = B[1,0] + B[1,1]*x + B[1,2]*x^2
    
    res = f.resultant(g)
    if res == n:
        listpoly.append([f,g])
    

    return listpoly



###########################################################
#twoquadratics                                            #
###########################################################
#Input:                                                   #
#  - two integers n and k                                 #
#                                                         #
#Output:                                                  #
#  - k pairs of quadratic polynomials                     #
#For each pair (f,g), f and g have a common root modulo n #
###########################################################

def twoquadratics(n, k):
    list = []
    p = 2
    compteur = 0

    while compteur < k:
        while kronecker(n,p) != 1 :
            p = next_prime(p)

        li = TQponctuel(n, p)
        for i in range(len(li)):
            list.append(li[i])

        compteur += len(li)
        p = next_prime(p)

    return list



##############################################################
#L2norm                                                      #
##############################################################
#Input:                                                      #
#  - two quadratic polynomials f and g                       #
#                                                            #
#Output:                                                     #
#  - the minimal value of the parameter integral             #
#    of V(s) = (F(x*s,y/s)*G(x*s,y/s))^2 over [-1;1]*[-1;1]  #
#  - the value s0 which minimizes T(s)                       #
##############################################################

def L2norm(f,g):

    a0 = f[0]
    a1 = f[1]
    a2 = f[2]
    b0 = g[0]
    b1 = g[1]
    b2 = g[2]

    #T(x) = V(x^(1/4))
    T(x) = 4/1575*(175*a2^2*b2^2*x^4 + 75*(a2^2*b1^2 + (2*a0*a2 + a1^2)*b2^2 + 2*(2*a1*a2*b1 + a2^2*b0)*b2)*x^3 + 175*a0^2*b0^2 + 63*(a0^2*b2^2 + 4*a1*a2*b0*b1 + a2^2*b0^2 + (2*a0*a2 + a1^2)*b1^2 + 2*(2*a0*a1*b1 + (2*a0*a2 + a1^2)*b0)*b2)*x^2 + 75*(2*a0^2*b0*b2 + a0^2*b1^2 + 4*a0*a1*b0*b1 + (2*a0*a2 + a1^2)*b0^2)*x)/x^2

    #DT(x) is the derivative of T(x). By finding its zero (which we assume to be unique), we can find the minimum of T(x)
    DT(x) = 4/1575*(700*a2^2*b2^2*x^3 + 150*a0^2*b0*b2 + 75*a0^2*b1^2 + 300*a0*a1*b0*b1 + 75*(2*a0*a2 + a1^2)*b0^2 + 225*(a2^2*b1^2 + (2*a0*a2 + a1^2)*b2^2 + 2*(2*a1*a2*b1 + a2^2*b0)*b2)*x^2 + 126*(a0^2*b2^2 + 4*a1*a2*b0*b1 + a2^2*b0^2 + (2*a0*a2 + a1^2)*b1^2 + 2*(2*a0*a1*b1 + (2*a0*a2 + a1^2)*b0)*b2)*x)/x^2 - 8/1575*(175*a2^2*b2^2*x^4 + 75*(a2^2*b1^2 + (2*a0*a2 + a1^2)*b2^2 + 2*(2*a1*a2*b1 + a2^2*b0)*b2)*x^3 + 175*a0^2*b0^2 + 63*(a0^2*b2^2 + 4*a1*a2*b0*b1 + a2^2*b0^2 + (2*a0*a2 + a1^2)*b1^2 + 2*(2*a0*a1*b1 + (2*a0*a2 + a1^2)*b0)*b2)*x^2 + 75*(2*a0^2*b0*b2 + a0^2*b1^2 + 4*a0*a1*b0*b1 + (2*a0*a2 + a1^2)*b0^2)*x)/x^3

    #We approximate the value of the zero of DT by dichotomy
    u1 = 1
    while (DT(u1) > 0):
        u1 = u1/2
    u2 = 1
    while (DT(u2) < 0):
        u2 = 2*u2
    while (u2 - u1 > 0.00000001):
        u3 = (u1+u2)/2
        if DT(u3) > 0:
            u2 = u3
        else:
            u1 = u3

    #s0 is the value that minimizes V(s)
    s0 = float(u1^(1/4))

    #norm is the minimum of V(s) = T(s^4)
    norm = float(T(u1))

    return [s0, float(norm)]


##############################################################
#rating                                                     #
##############################################################
#Input:                                                      #
#  - two quadratic polynomials f and g                       #
#                                                            #
#Output:                                                     #
#  - alpha(f, 2000)                                          #
#  - alpha(g, 2000)                                          #
#  - the rating of the pair (f,g)                            #
#  - L2norm(f,g)                                             #
##############################################################

def rating(f, g):
    a = alpha(f, 2000)
    b = alpha(g, 2000)
    return [a, b, a + b + 0.5*log(L2norm(f,g)[1]), L2norm(f,g)[0]]



#####################################################################
#PolynomialSelection                                                #
#####################################################################
#Input:                                                             #
#  - two integers n and k                                           #
#                                                                   #
#Output:                                                            #
#  - a pair of polynomials (f,g), along with its rating             #
#    this pair is chosen among k pairs for having the lowest rating #       
#####################################################################

def PolynomialSelection(n, k):
    Liste = twoquadratics(n,k)
    min = float(100)

    notabis = []

    for i in range(k):
        print ""
        print i
        nota = rating(Liste[i][0], Liste[i][1])
        print Liste[i]
        print nota
        if nota[2] < min:
            print "on a un nouveau minimum"
            f = Liste[i][0]
            g = Liste[i][1]
            min = nota[2]
            notabis = nota
    print ""
    print "le choix optimal est "
    print "f = ", f
    print "g = ", g
    print "alpha(f, 2000) = ", notabis[0]
    print "alpha(g, 2000) = ", notabis[1]
    print "notation(f,g) = ", notabis[2]
    print "s optimal = ", notabis[3]
