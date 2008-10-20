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
    assert(u0==0)
    assert(v0==0)
    space=u1*v1
    print "Rotation space: %r" % space 
    rdict=dict(umax=u1,vmax=v1,f=f,g=g)
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
    df=f.derivative()
    dg=g.derivative()
    ddf=df.derivative()
    ddg=dg.derivative()
    rdict['fdg_gdf']=f*dg-g*df
    return rdict

def rotation_clear(rdict):
    sarr=rdict['sarr']
    for i in range(len(sarr)): sarr[i]=0.0

def extend_ppow_list(rdict,l):
    ppow = rdict['ppow']
    ppows = len(ppow)
    pmaxi = ppow[ppows-1]
    while ppows < l+1:
        pmaxi *= ppow[1]
        ppow.append(pmaxi)
        ppows+=1


def rotation_handle_p(rdict,p):
    u0=v0=0
    ff,gg=rdict['f'],rdict['g']
    f,g=rdict['f'],rdict['g']
    x=Integers()['x'].gen()
    phi=x
    # The denominator exponent here is only controlled by the number of
    # fixed coeffs of l. At the beginning, l has zero fixed coeffs, for
    # which scale0 holds.
    scale0=float(log(p))*p/(p+1)
    l0,ld,md,m,k=0,0,0,0,0
    scale = scale0
    ZP3.<l,u,v>=ZZ['l','u','v']
    rdict['h']=f(l)+(u*l+v)*g(l)
    rdict['history']=[]
    rdict['logbook']=[]
    rotation_inner(rdict,p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,0)

def rotation_inner(rdict,p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,d):
    maxpow = 20
    if m >= maxpow:
        return

    rdict['history'].append((p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,d))
    # p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,d=rdict['history'][-1]

    callnumber=len(rdict['history'])-1
    prefix = " " * d
    msg=prefix+"Entering recursive call number %d, depth %d" % (callnumber,d)
    rdict['logbook'].append(msg)
    print msg

    a = vector([p^k,0])
    b = vector([0,p^k])
    uv0=vector((u0,v0))

    K=GF(p); KP=K['x']; KR=FractionField(KP)
    Z=Integers(); ZP=Z['x'];
    x=ZP.gen()
    sarr=rdict['sarr']
    fence=vector([rdict['umax'],rdict['vmax']])

    h = rdict['h']
    l,u,v=h.parent().gens()
    nh = h(phi(l),u0+p^k*u,v0+p^k*v)
    assert valuation(nh.content(),p)==m
    assert nh/p^m == ff(l)+(u*phi(l)+v)*gg(l)
    assert valuation(nh.derivative(l).content(),p)>=m
    # assert valuation(nh.derivative(l).content(),p)==ld+md

    # Our expression is c0(l)+u*cu(l)+v*cv(l)
    c0 = ZP([c % p^maxpow for c in ff.coeffs()])
    cv = ZP([c % p^maxpow for c in gg.coeffs()])
    cu = phi * cv
    assert valuation((nh - p^m*(c0(l)+u*cu(l)+v*cv(l))).content(),p) >= maxpow
    assert cv == Integers(p^maxpow)['x'](rdict['g'](phi))

    # The assertions on the trivariate polynomial nh above contain more
    # that this, but here follows a breakdown just in case.
    assert p^m*ff(l) == rdict['f'](phi(l))+((u0*phi(l)+v0)*rdict['g'](phi(l)))
    assert p^m*(u*phi(l)+v)*gg(l)==p^k*((u*phi(l)+v)*rdict['g'](phi(l)))
    assert p^m*gg(l)==p^k*rdict['g'](phi(l))

    cc0 = valuation(c0.content(),p)
    ccu = valuation(cu.content(),p)
    ccv = valuation(cv.content(),p)
    assert ccu >= ccv

    c = min(cc0,ccu,ccv)
    assert c == 0

    # zview(mround(matrix(sbound,sbound,sarr)-complete),u0,v0,p^k)

    dc0 = c0.derivative()
    dcu = cu.derivative()
    dcv = cv.derivative()

    # Now we know that the minimum valuation is zero.

    g0 = KP(c0)
    gu = KP(cu)
    gv = KP(cv)

    if not g0.is_constant():
        # Then it's the typical case, as encountered for instance at the
        # beginning of the root sieve.
        # Each possible l value must be tried, and will give rise to
        # potentially intersecting lattices. 
        for l in range(p):
            if k==0: l0=l
            g0l = g0(l); dg0l = KP(dc0)(l)
            gul = gu(l); dgul = KP(dcu)(l)
            gvl = gv(l); dgvl = KP(dcv)(l)
            if gvl == 0: continue
            msg=prefix+"Looking at roots phi(%r), where phi==%r" %(l,phi)
            rdict['logbook'].append(msg)
            print msg
            igl = 1/gvl
            du=0
            dv0=Z(-g0l * igl)
            dv=dv0
            u1 = u0 + p^k * du
            v1 = v0 + p^k * dv
            uv1=vector((u1,v1))
            rdict['logbook'].append((uv1, a-phi(l)*b, p* b, -scale/(p-1)))
            hits = light_array(sarr, uv1, a-phi(l)*b, p* b, fence,-scale/(p-1))
            msg=prefix+ "%d hits" %hits
            rdict['logbook'].append(msg)
            print msg
            if hits == 0: continue
            # We're optimistic here ; we're counting on the fact that we
            # expect the derivative not to cancel.
            
            # However, there is the possibility that we reach the
            # cancellation point. This only once per (p,p) square,
            # meaning that u,v are constrained.

            msg=prefix+ "Looking for multiple roots at phi(%r), where phi==%r" %(l,phi)
            rdict['logbook'].append(msg)
            print msg
            # The good equation for u is given by:
            df=rdict['f'].derivative()
            dg=rdict['g'].derivative()

            # Try to guess the right value for du.
            Zm1=Integers(p^(m+1))
            # rhs = Z(Zm1((f*dg-df*g)(phi(l))/cv(l)^2-u0))
            rhs = Z(K(((f*dg-g*df)(phi(l))/cv(l)^2-u0)/p^(m-ld)))
            # va = valuation(rhs,p)
            # assert va >= md
            # assert k >= md
            # if va < k:
                # # Then by no means the derivative will vanish.
                # continue
            range_du=[]
            if k-m+ld > 0 and rhs == 0:
                # Then all u's will do.
                range_du=range(p)
            else:
                range_du = [rhs]
            ## u = Z(-igl*K((((df+(u0*x+v0)*dg+u0*g)(phi))(l)-ff(l)*(p^k*dg(phi))(l))/p^m))
            ## du = Z((-dg0l * gvl + g0l * dgvl) * igl ^ 2)
            ## u = Z(Z(Integers(p^(m+1))((f*dg-df*g)(phi(l))/cv(l)^2)-uv0[0])/p^k)
            for du in range_du:
                dv = (dv0 - du *l0) % p
                u1 = u0 + p^k * du
                v1 = v0 + p^k * dv
                uv1=vector((u1,v1))
                rdict['logbook'].append((uv1, p*a, p*b, scale/(p-1)/p))
                hits = light_array(sarr,
                            uv1, p*a, p*b, fence, scale/(p-1)/p)
                msg=prefix+ "%d hits" %hits
                rdict['logbook'].append(msg)
                print msg
                if hits == 0:
                    continue
                # scale/(p-1)/p is scale/(p-1)-scale/p ; hence it brings
                # back the cells to the contribution -scale/p, which is the
                # contribution from a _single_ zero at this location, not
                # liftable because of ramification. Later recursive calls
                # will investigate the possibility that despite the multiple
                # root mod p, we still get roots at higher precision.
                nphi=phi(l+p*x)
                nff=ZP((ff(l+p*x)+(du*nphi+dv)*gg(l+p*x))/p)
                ngg=gg(l+p*x)
                ## l,u,v=h.parent().gens()
                ## nnh=nh(nphi(l),uv1[0]+p*u,uv1[1]+p*v)
                ## assert nnh == p*(nff(l)+(u*nphi(l)+v)*ngg(l))
                # ff,gg,uv0,l0,m,k,scale,phi=nff,ngg,uv1,l0,m+1,k+1,scale/p,nphi
                rotation_inner(rdict,p,nff,ngg,u1,v1,l0,ld+1,md+1,m+1,k+1,scale/p,nphi,d+1)
    elif ccv > 0:
        # ccv is zero, so ccu is even larger, which implies that cc0 has
        # to be zero. Since c0 is a constant polynomial, this implies
        # that there is no solution. So we finish here.
        return
    else:
        msg=prefix+ "We are here in recursive call number %d" % callnumber
        rdict['logbook'].append(msg)
        print msg
        # cc0 is not zero, which means that our expression has an extra
        # root only for specific u,v pairs. Basically, we have a linear
        # term in u and v which must be cancelled.
        constant_term=g0.constant_coefficient()
        dv0=Z(-constant_term/gv(l0))
        dv=dv0
        u1 = u0
        v1 = v0 + p^k * dv
        uv1=vector((u1,v1))
        assert ccv==0
        # Furthermore, our expression is c0(l)+u*cu(l)+v*cv(l) ; We know
        # that cu == phi * cv, hence our equation is simply v == -phi * u
        msg=prefix+ "Looking for places where the multiple root at %r may lift"%phi
        rdict['logbook'].append(msg)
        print msg
        rdict['logbook'].append((uv1, a-l0*b, p*b,-scale))
        hits=light_array(sarr, uv1, a-l0*b, p* b, fence,-scale)
        msg=prefix+ "%d hits" %hits
        rdict['logbook'].append(msg)
        print msg
        if hits == 0:
            return
        msg=prefix+"%r"%(K['l','u','v'](c0(l)+u*cu(l)+v*cv(l)))
        rdict['logbook'].append(msg)
        print msg
        for du in range(p):
            dv = (dv0-du*l0) % p
            uvscal=vector((du,dv))
            msg=prefix+ "Looking for lift specifically at u,v === %r" % uvscal
            rdict['logbook'].append(msg)
            print msg
            u2 = u0 + p^k * du
            v2 = v0 + p^k * dv
            uv2=vector((u2,v2))
            msg=prefix+ "global sub-lattice %r+%r,%r" % (uv2,p*a,p*b)
            rdict['logbook'].append(msg)
            print msg
            # We know that phi === l0 mod p, so (u0*phi+v0) is zero mod
            # p. ff itself is also zero mod p. However, some subtleties
            # can occur in the difference between the two.
            nff=ZP((ff+(du*phi+dv)*gg)/p)
            l,u,v=h.parent().gens()
            assert nh(l,du+p*u,dv+p*v)/p^(m+1) == nff(l)+(u*phi(l)+v)*gg(l)
            ## ff,u0,v0,m,k=nff,u1,v1,m+1,k+1
            rotation_inner(rdict,p,nff,gg,u2,v2,l0,ld,md,m+1,k+1,scale,phi,d+1)
    msg=prefix+ "Finishing recursive call number %d" % callnumber
    rdict['logbook'].append(msg)
    print msg


def filter_logbook(rdict,u,v):
    vec=vector((u,v))
    for s in rdict['logbook']:
        if type(s) == type(''):
            print s
        else:
            st="sub-lattice %r+%r,%r: %f" % s
            if ((vec-s[0])*matrix([s[1],s[2]])^-1) in vec.parent():
                st="** " + st
            print st

def filter_logbook_quick(rdict,u,v):
    vec=vector((u,v))
    for s in rdict['logbook']:
        if type(s) != type(''):
            if ((vec-s[0])*matrix([s[1],s[2]])^-1) in vec.parent():
                st="sub-lattice %r+%r,%r: %f" % s
                st="** " + st
                print st


def light_array(sarr,x0,a,b,fence,value):
    x = copy(x0)
    pos = x[0] * fence[1]
    assert b[0] == 0
    dpos = a[0] * fence[1]
    hits=0
    while x[0] < fence[0]:
        # reduce within the row, since we know that b[0] is zero
        x[1] = x[1] % b[1]
        while x[1] < fence[1]:
            sarr[pos + x[1]] += value
            hits+=1
            x[1] += b[1]
        x += a
        pos += dpos
    return hits

# These are helpers, to be removed once the code works ok.
def zview(m,i0,j0,d):
    ncols = ((m.ncols() - j0) / d).ceil()
    nrows = ((m.nrows() - i0) / d).ceil()
    n=matrix(nrows,ncols)
    for i in range(0,nrows):
        for j in range(0,ncols):
            n[i,j]=m[i0+i*d,j0+j*d]
    return n

def mround(m):
    d=matrix(ZZ,m.nrows(),m.ncols())
    for k in range(m.nrows()):
        for l in range(m.nrows()):
            d[k,l]=ZZ(floor(m[k,l]*1000))
    return d

def mprint(m):
    print m.str()
