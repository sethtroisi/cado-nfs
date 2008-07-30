attach "alpha.sage"
attach "rotation_bound.sage"

# for reading kleinjung.res
import re

def lemme_21(n,d,ad,p,m):
    """Implements Kleinjung's lemma 2.1. ad must satisfy ad*m^d==n mod p.
    ad may also be equal to zero, in which case it is recovered
    automatically"""
    if ad==0:
        ad=ZZ(Integers(p)(n/m^d))
    assert(ad*Integers(p)(m)^d-n==0)
    coeffs=[ad]
    deltas=[]
    rs=[]
    r=n
    a=ad
    Z=Integers()
    mm=[1]
    immp=[Integers(p)(1)]
    imp=1/Integers(p)(m)
    for i in [1..d]:
        mm.append(mm[-1]*m)
        immp.append(immp[-1]*imp)
    for i in reversed([0..d-1]):
        r=Z((r-a*mm[i+1])/p)
        ai_mod_p=Z(r*immp[i])
        k=(r/mm[i]-ai_mod_p)/p
        kr=round(k)
        delta=p*(kr-k)
        a=p*kr+ai_mod_p
        rs.append(r)
        coeffs.append(a)
        deltas.append(delta)
    coeffs.reverse()
    rs.reverse()
    deltas.reverse()
    f=Integers()['x'](coeffs)
    g=Integers()['x']([-m,p])
    #return (f,deltas,rs)
    return (f,g)

# For example:

# n=1230186684530117755130494958384962720772853569595334792197322452151726400507263657518745202199786469389956474942774063845925192557326303453731548268507917026122142913461670429214311602221240479274737794080665351419597459856902143413
# f=ZZ['x']([-277565266791543881995216199713801103343120,-18185779352088594356726018862434803054,6525437261935989397109667371894785,-46477854471727854271772677450,-5006815697800138351796828,1276509360768321888,265482057982680])
# g=ZZ['x']([-1291187456580021223163547791574810881,34661003550492501851445829])
# f0,_=lemme_21(n,6,f[6],g[1],-g[0])
# (f-f0)/g

# Do for instance:
# c136=3835722565249558322197842586047190634357074639921908543369929770877061053472916307274359341359614012333625966961056935597194003176666977
# n=c136
# d=5 
# tuples=import_kleinjung_dot_res(n,d,10^19,"/net/tiramisu/cado1/cado/Examples/35_143+_c136/kleinjung.res")

def translate_polypair_string(n,d,normbound,s):
    """Given a report string of the form ``p=... m=... norm=...'', return
    the corresponding polynomial pair if norm<normbound, or None
    otherwise"""
    ma=re.match("^p=(\d+) m=(\d+)(?: norm=([\d\.]+e\+\d+))?",s)
    # No match is something normal
    if ma == None: return
    tu=ma.groups()
    f=None
    g=None
    if tu[2] == None:
        p=ZZ(tu[0])
        m=ZZ(tu[1])
        f,g=lemme_21(n,d,0,p,m)
        norm=best_l2norm_tk(f)
        if norm >= normbound: return
    else:
        norm=RR(tu[2])
        if norm >= normbound: return
        p=ZZ(tu[0])
        m=ZZ(tu[1])
        f,g=lemme_21(n,d,0,p,m)
    s0="p=%r m=%r supnorm=%.2e"%(p,m,norm)
    print s0
    return (f,g)

def import_kleinjung_dot_res(n,d,normbound,filename):
    """This function reads and parses the kleinjung.res file. Only
    reports below the specified normbound are processed."""
    fi=open(filename,"r")
    tuples=[translate_polypair_string(n,d,normbound,s) for s in fi.readlines()]
    tuples=[t for t in tuples if t != None]
    fi.close()
    return tuples

def optimize_tuple(f,g):
    """Optimizes the polynomial pair f,g. Return the best logmu+alpha"""
    p=g[1]
    m=-g[0]
    s=skew_l2norm_tk(f)
    a0=alpha_projective(f,2000)
    a1=alpha_affine(f,2000)
    y0=flog(best_l2norm_tk(f))+a0+a1
    toto=lognorm_plus_alpha_rot_scons_linear
    c=[(a0+toto(f,g,l2norm_tk,s,i),i) for i in [0..14]]
    cm=min(c)
    s0= "p=%r m=%r"%(p,m)
    s1= "%.2f -> %.2f [10^%d]" % (y0,cm[0],cm[1])
    print s0 + " " + s1
    return cm

def optimize_tuples(l):
    """Given a list as returned by import_kleinjung_dot_res, optimizes
    all one by one"""
    cms=[(optimize_tuple(l[i][0],l[i][1]),i) for i in range(len(l))]
    mi=min(cms)
    i=mi[1]
    f=l[i][0]
    g=l[i][1]
    y=mi[0][0]
    w=mi[0][1]
    print "min: [#%d] p=%r m=%r %.2f [10^%d]" % (i,g[1],-g[0],y,w)

def rotation (n,d,p,m):
    """Computes optimal rotation for degree-d polynomial for integer n with
       linear polynomial p*x-m"""
    B = 2000
    f,g = lemme_21(n,d,0,p,m)
    x = f.parent().gen()
    s = skew_l2norm_tk (f)
    w = 2
    r = lognorm_plus_alpha_rot_scons_linear (f, g, l2norm_tk, s, flog10(w))
    oldr = Infinity
    while r < oldr:
        oldr = r
        w = 2 * w
        r = lognorm_plus_alpha_rot_scons_linear (f, g, l2norm_tk, s, flog10(w))
    # the constant polynomial varies up to V = sqrt(w*s), and the linear
    # polynomial up to U = sqrt(w/s)
    U = ceil(sqrt(w/s))
    V = ceil(sqrt(w*s))
    print "rotation bounds:", U, V
    sys.stdout.flush()
    rotation_inner(f,g,range(-U,U+1),range(-V,V+1))

def rotation_inner(f,g,u_range, v_range):
    """returns the multiplier lambda=u*x+v giving the best yield for
    f+lambda*g ; u and v are prescribed to the ranges u_range and
    v_range"""
    B=2000
    x = f.parent().gen()
    lognorm = flog (best_l2norm_tk(f))
    a0=alpha_projective(f,B)
    a1=alpha_affine(f,B)
    a=a0+a1
    E = lognorm + a
    print "original polynomial: lognorm=%.2f alpha=%.2f E=%.2f" % (lognorm,a,E)
    Emin = E
    umin = 0
    vmin = 0
    suma1=0
    suma1_squares=0
    for u in u_range:
        for v in v_range:
            frot = f + (u*x+v)*g
            lognorm = flog (best_l2norm_tk (frot))
            a1 = alpha_affine (frot, B)
            a = a0 + a1
            E = lognorm + a
            suma1+=a1
            suma1_squares+=a1*a1
            if E < Emin:
                Emin = E
                umin = u
                vmin = v
                print "u=%r v=%r lognorm=%.2f alpha=%.2f E=%.2f" % (u,v,lognorm,a,E)
                sys.stdout.flush()
    me=suma1/len(u_range)/len(v_range)
    sd=sqrt(suma1_squares/len(u_range)/len(v_range)-me^2)
    print "Observed alpha_affine: mean=%.2f, sdev=%.2f" % (me,sd)
    return p, m, umin*x+vmin, Emin





def sylvestermatrix(p1,p2):
    m=matrix(p1(0).parent(),p1.degree()+p2.degree())
    for i in range(p1.degree()):
        for j in [0..p2.degree()]:
            m[i,i+j]=p2[j]
    for j in range(p2.degree()):
        for i in [0..p1.degree()]:
            m[j+p1.degree(),i+j]=p1[i]
    return m

def discriminant_polynomial(f,g):
    t0=cputime()
    # We want the discrimnant in x of (f+(ux+v)g). Sage can't do that, so
    # we'll resort on computing the Sylvester matrix ourselves.
    ZP3.<x,u,v>=Integers()['x','u','v']
    frot=(f(x)+(u*x+v)*g(x)).polynomial(x)
    m=sylvestermatrix(frot, frot.derivative())
    sys.stdout.write("Computing discriminant polynomial...");
    sys.stdout.flush()
    t0=cputime()
    disc=m.determinant()
    z=disc.monomials()
    lc=frot[frot.degree()]
    d1=sum([zz*ZZ(QQ(disc.monomial_coefficient(zz)/lc)) for zz in z])
    t1=cputime()
    sys.stdout.write("done in %.2fs\n" % (t1-t0))
    sys.stdout.flush()
    return d1

def discriminant_polynomial_pari(f,g):
    """This is instantaneous. However, once we've done this, there's no
    coming back"""
    ZP3=PolynomialRing(ZZ,['x','u','v'])
    x,u,v=ZP3.gens()
    dd=pari(f(x)+(u*x+v)*g(x)).poldisc(x)
    return dd

def rotation_global_contrib_commonroots(rdict):
    """Fills in the global_contrib_commonroots entry in the provided
    rotation dictionary. This makes sense only for SNFS"""
    global_contrib=float(0.0)
    f=rdict['f']
    g=rdict['g']
    B=rdict['B']
    res=f.resultant(g)
    while true:
        p=trial_division(res,B)
        if p==res:
            break
        scale=p*float(log(p))/(p^2-1)
        global_contrib-=scale
    for i in range(len(rdict['sarr'])):
        rdict['sarr'][i]+=global_contrib
    rdict['global_contrib_commonroots']=global_contrib

def rotation_init(f,g,u0,u1,v0,v1):
    """Sets up rotation for using the area [u0..u1] x [v0..v1]"""
    lu=u1-u0+1
    lv=v1-v0+1
    space=lu*lv
    print "Rotation space: %r" % space 
    rdict=dict(u0=u0,u1=u1,v0=v0,v1=v1,f=f,g=g)
    rdict['sarr']=[float(0.0) for i in range(space)]
    # Estimate the common value for the different cells.
    B=2000
    rdict['B']=B
    go=sum([float(log(p))/(p-1) for p in prime_range(B)])
    # for i in range(len(rdict['sarr'])): rdict['sarr'][i]+=go
    rdict['global_offset']=go
    # Also the contribution from common roots (SNFS only)
    rotation_global_contrib_commonroots(rdict)
    gcc=rdict['global_contrib_commonroots']
    rdict['fgquo_prime']=((-f/g).derivative())
    # This gives the exceptional values for u (once divided by g(l)^2)
    rdict['fggf']=f*g.derivative()-g*f.derivative()
    return rdict

def lift_root(f,l,p,k):
    x=f.parent().gen()
    fl=f(l)
    v=valuation(fl,p)
    assert v>=k
    g=f.parent()(f(l+p^k*x)/p^v)
    R=GF(p)
    RP=R['x']
    gg=RP(g)
    if gg.degree() != 1:
        print "Cannot lift further ; degree == %d" %gg.degree()
        return
    e=-gg[0]/gg[1]
    print "Valuation %r^%r: %r" %(p,k,e)
    if k < 10:
        l+=p^k*e
        lift_root(f,l,p,k+1)

def rotation_handle_p(rdict,p):
    scale=float(log(p))/(p+1)
    rotation_handle_p_recursive(rdict,p,0,0,0,0,float(log(p))/(p+1))

def rotation_handle_p_recursive_inner(rdict,p,l0,u0,v0,k):
    """Recursively computes alpha contributions mod p for (l,u,v)
    congruent to (l0,u0,v0) mod p^k"""
    # If we denote S(l,u,v)=(f(l)+(u*l+v)*g(l)), then Sl,Su,Sv the partial
    # derivatives wrt l,u,v, and finally X0 the triple (l0,u0,v0), then
    # we have forced:
    # S(X0)=0 mod p^k
    # Sl(X0)=0 mod p^k -- because otherwise we're talking simple roots,
    # which are easily worked out.
    #
    # The question is how does S(X0+p^kX1) behave, for X=(l1,u1,v1)) ?
    # We have:
    # S(X0+p^kX1)=S(X0)+p^k(l1*Sl(X0)+u1*Su(X0)+v1*Sv(X0))+O(p^(2k))
    #            =p^k*(S(X0)/p^k+u1*Su(X0)+v1*Sv(X0)) + O(p^(2k))
    #
    # This implies that our congruences can be refined considerably. For
    # each value u1 in [0..p^k-1], we have a corresponding value v1 such
    # that the equation above yields cancellation modulo p^(2k). Note
    # that:
    #
    # - we want to count contributions also modulo p^(k+1) and so on.
    # This will give more values of v1, many of which will give
    # ``terminating'' congruences. This is done in k different passes
    # through the v-line.
    #
    # - Su(X0) and Sv(X0) are in fact pretty similar. Respectively
    # l0*g(l0) and g(l0). g(l0) is assumed to be a unit mod p, so it has
    # a computable inverse modulo p^(2k). This yields the equation:
    #
    # v1 = -u1 * l0 - (1/g(l0)) * (S(X0)/p^k)
    #
    # as in the first-order case, forcing the derivative to vanish modulo
    # p^k yields specific values of u, since:
    # Sl = f'(l)+(ul+v)g'(l)+ug(l)
    K=GF(p)
    pl=[0 for i in [0..k]]
    pl[0]=p^k
    for i in [1..k]:
        pl[i]=pl[i-1]*p

    pk=pl[0]
    p2k=pl[k]
    R0=Integers(pk)
    R=Integers(p2k)
    RP.<x>=R['x']
    fl0=RP(rdict['f'])(l0)
    gl0=RP(rdict['g'])(l0)
    assert(gl0.is_unit())
    # igl0 is needed to high precision in order to compute the
    # exceptional u1 value. For the congruence equation, it's unused.
    igl0=1/R(gl0)
    SX0=fl0+(u0*l0+v0)*gl0
    # How does the congruence equation look like ?
    veq_0=-R0(igl0)*R0(Integers()(Integers()(SX0)/pk))
    veq_1=R0(-l0)
    # So basically it's v1 = veq_0 + veq_1 * u1

    # The global u range is rdict.u0..rdict.u1
    # We want u0+pk*u1 in this range.
    u1min = ceil((rdict['u0']-u0)/pk)
    u1max = floor((rdict['u1']-u0)/pk)
    v1min = ceil((rdict['v0']-v0)/pk)
    v1max = floor((rdict['v1']-v0)/pk)

    u1 = u1min

    lv=rdict['v1']-rdict['v0']+1
    lvpk=lv*pk
    v1p = veq_0 + veq_1 * u1

    uoffset = (u0+u1min*pk - rdict['u0']) * lv
    voffset = (v0+v1min*pk - rdict['v0'])

    u1_exceptional = R0(Integers()(RP(rdict['fggf'])(l0) * igl0^2 - u0)/p)

    u1p = R0(u1)

    # Pre-compute v1 mod p^l for l in [1..k] ; we leave cell [0] unused.
    v1ls=[0 for i in [0..k]]
    v1ls[k]=Integers()(v1p)
    for l in reversed([1..k-1]):
        v1ls[l]=v1ls[l+1] % pl[l]

    while u1 <= u1max:
        # Hit points in [vmin,vmax] that correspond to v1...

        # For congruences to powers less than p^(2k-1), we always have
        # contributions. We sum them up in the most pedestrian way. A
        # priori, congruences will stop. So the contribution at each
        # prime power is only 1/p^(k+l)

        for l in [1..k-1]:
            # modulo p^(k+l)
            # find v in [v1min..v1max] such that v==v1 mod p^l
            dv=(v1ls[l]-v1min)%pl[l]
            z=offset+voffset+pk*dv
            dv+=v1min
            while dv < v1max:
                rdict['sarr'][z]-=1/pl[l]
                z+=pl[l]
                dv+=pl[l]

        if u1p != u1_exceptional:
            # The generic case. Congruences modulo p^(2k) are
            # liftable to arbitrary precision since the derivative does
            # not vanish. Hence the contribution to the valuation is
            # 1/(p^2k-1).
            l=k
            dv=(v1ls[l]-v1min)%pl[l]
            z=offset+voffset+pk*dv
            dv+=v1min
            while dv < v1max:
                rdict['sarr'][z]-=1/(pl[l]-l)
                z+=pl[l]
                dv+=pl[l]
        else:
            # Then we need to recurse, yeah !
            # We do count the p^2k contribution right now. Higher order
            # contributions will be handled in the recursive calls.
            # contributions are 
            l=k
            dv=(v1ls[l]-v1min)%pl[l]
            z=offset+voffset+pk*dv
            dv+=v1min
            while dv < v1max:
                rdict['sarr'][z]-=1/pl[l]
                rotation_handle_p_recursive_outer(rdict,p,l0,u0,v0,2*k)
                z+=pl[l]
                dv+=pl[l]

        u1+=1
        u1p+=1
        v1ls[k]=Integers(v1ls[k]+veq_1)
        for l in reversed([1..k-1]):
            v1ls[l]=v1ls[l+1] % pl[l]
        offset+=lv*pk


def rotation_handle_p_recursive_outer(rdict,p,l0,u0,v0,k):
    """l0 is known to be a root at (u0,v0) mod p^k"""
    ### This code is not finished.
    assert(false)
    K=GF(p)
    KP=K['x']
    x=KP.gen()
    KR=FractionField(KP)
    scale0=float(log(p))/(p-1)
    scale=p*float(log(p))/(p^2-1)
    fgquo_p=KR(rdict['fgquo_prime'])
    fp=KP(rdict['f'])
    gp=KP(rdict['g'])
    sarr=rdict['sarr']
    for l in K:
        fl=fp(l)
        gl=gp(l)
        # Then we want v == -fl/gl - u*l mod p if ever gl == 0, then we
        # want fl==0 and that's all (and independent of u and v) ; if
        # this is the case, it would mean that f and g have a common root
        # mod p. It is handled as a global contribution further up (and
        # only for SNFS)
        if gl==0: continue
        zl=-fl/gl
        # print "l=%r -> z=%r" %(l,zl)
        # For which u values are multiple roots possible ?
        u_excep=fgquo_p(l)
        # Update array.
        pos=0
        v0=rdict['v0']
        v1=rdict['v1']
        u=rdict['u0']
        u1=rdict['u1']
        up=K(u)
        lv=v1-v0+1
        while u <= u1:
            # want v0+k == zl-ul mod p
            k=Integers()(K(zl-u*l-v0))
            if up != u_excep:
                print "root at l=%r for u,v = %r,%r" %(l,u,v0+k)
                while k < lv:
                    sarr[pos+k]-=scale
                    k+=p
            else:
                # We know that l will be a multiple root.
                print "Multiple root at l=%r for u,v = %r,%r" %(l,u,v0+k)
                while k < lv:
                    # Count only the local contribution for now. Higher
                    # contributions will be counted separately.
                    sarr[pos+k]-=scale0
                    k+=p
            pos+=lv
            u+=1
            up+=1

def rotation_handle_p(rdict,p):
    K=GF(p)
    KP=K['x']
    x=KP.gen()
    KR=FractionField(KP)
    scale0=float(log(p))/(p-1)
    scale=p*float(log(p))/(p^2-1)
    fgquo_p=KR(rdict['fgquo_prime'])
    fp=KP(rdict['f'])
    gp=KP(rdict['g'])
    sarr=rdict['sarr']
    for l in K:
        fl=fp(l)
        gl=gp(l)
        # Then we want v == -fl/gl - u*l mod p if ever gl == 0, then we
        # want fl==0 and that's all (and independent of u and v) ; if
        # this is the case, it would mean that f and g have a common root
        # mod p. It is handled as a global contribution further up (and
        # only for SNFS)
        if gl==0: continue
        zl=-fl/gl
        # print "l=%r -> z=%r" %(l,zl)
        # For which u values are multiple roots possible ?
        u_excep=fgquo_p(l)
        # Update array.
        pos=0
        v0=rdict['v0']
        v1=rdict['v1']
        u=rdict['u0']
        u1=rdict['u1']
        up=K(u)
        lv=v1-v0+1
        while u <= u1:
            # want v0+k == zl-ul mod p
            k=Integers()(K(zl-u*l-v0))
            if up != u_excep:
                print "root at l=%r for u,v = %r,%r" %(l,u,v0+k)
                while k < lv:
                    sarr[pos+k]-=scale
                    k+=p
            else:
                # We know that l will be a multiple root.
                print "Multiple root at l=%r for u,v = %r,%r" %(l,u,v0+k)
                while k < lv:
                    # Count only the local contribution for now. Higher
                    # contributions will be counted separately.
                    sarr[pos+k]-=scale0
                    k+=p
            pos+=lv
            u+=1
            up+=1

#def rotation_using_sieve2(f,g,u_range, v_range,disc):
#    #sarr=rotation_using_sieve(f,g,u_range, v_range)
#    B=2000
#    fraction=0
#    for p in prime_range(B):
#        KP2=GF(p)['u','v']
#        KP1.<v>=GF(p)['v']
#        dp=KP2(disc)
#        # This will give us the update points on one of the sieve lines
#        for u in u_range:
#            dpu=KP1(dp(u,v))
#            if dpu==0:
#                # print "disc(%d,x)==0 mod %d"%(u,p)
#                fraction+=flog(p)^2
#                continue
#            rdisc=dpu.roots()
#            # print "%d roots of the discriminant mod %d" % (len(rdisc),p)
#            fraction+=len(rdisc)*float(p)^2/p
#    return fraction
#    #return sarr
