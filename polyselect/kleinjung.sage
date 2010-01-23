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
    a=[0 for i in [0..d]]
    a[d] = ad
    r=[ 0 for i in [0..d] ]
    r[d]=n
    delta=[ 0 for i in [0..d] ]
    Z=Integers()
    mm=[1]
    immp=[Integers(p)(1)]
    imp=1/Integers(p)(m)
    for i in [1..d]:
        mm.append(mm[-1]*m)
        immp.append(immp[-1]*imp)
    for i in reversed([0..d-1]):
        r[i]=Z((r[i+1]-a[i+1]*mm[i+1])/p)
        ai_mod_p=Z(r[i]*immp[i])
        k=(r[i]/mm[i]-ai_mod_p)/p
        kr=Integers()(round(k))
        delta[i]=p*(kr-k)
        a[i]=p*kr+ai_mod_p
#        deltas.append(delta)
        assert n == sum([a[d-j]*m^(d-j)*p^j for j in range(0,d-i)])+r[i]*p^(d-i)
    if a[0] != r[0]:
        print "warning, using r0 instead of r0 (lemma error)"
        a[0]=r[0]
    f=Integers()['x'](a)
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
    #rotation_inner(f,g,range(-U,U+1),range(-V,V+1))

# def rotation_global_contrib_commonroots(rdict):
#     """Fills in the global_contrib_commonroots entry in the provided
#     rotation dictionary. This makes sense only for SNFS"""
#     global_contrib=float(0.0)
#     f=rdict['f']
#     g=rdict['g']
#     B=rdict['B']
#     res=f.resultant(g)
#     while true:
#         p=trial_division(res,B)
#         if p==res:
#             break
#         scale=p*float(log(p))/(p^2-1)
#         global_contrib-=scale
#     for i in range(len(rdict['sarr'])):
#         rdict['sarr'][i]+=global_contrib
#     rdict['global_contrib_commonroots']=global_contrib
# 
def rotation_init(f,g,u0,u1,v0,v1):
    """Sets up rotation for using the area [u0..u1] x [v0..v1]"""
    assert(u0==0)
    assert(v0==0)
    space=u1*v1
    #print "Rotation space: %r" % space 
    rdict=dict(umax=u1,vmax=v1,f=f,g=g)
    rdict['sarr']=[float(0.0) for i in range(space)]
    # Estimate the common value for the different cells.
    B=2000
    rdict['B']=B
    go=sum([float(log(p))/(p-1) for p in prime_range(B)])
    # for i in range(len(rdict['sarr'])): rdict['sarr'][i]+=go
    rdict['global_offset']=go
    # Also the contribution from common roots (SNFS only)
    ## rotation_global_contrib_commonroots(rdict)
    ## gcc=rdict['global_contrib_commonroots']
    rdict['common_contribution']=0
    rdict['fgquo_prime']=((-f/g).derivative())
    #print rdict['fgquo_prime']    
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


# def reference_rotation_handle_p(rdict,p):
#     u0=v0=0
#     ff,gg=rdict['f'],rdict['g']
#     f,g=rdict['f'],rdict['g']
#     x=Integers()['x'].gen()
#     phi=x
#     # The denominator exponent here is only controlled by the number of
#     # fixed coeffs of l. At the beginning, l has zero fixed coeffs, for
#     # which scale0 holds.
#     scale0=float(log(p))*p/(p+1)
#     l0,ld,md,m,k=0,0,0,0,0
#     scale = scale0
#     ZP3.<l,u,v>=ZZ['l','u','v']
#     rdict['h']=f(l)+(u*l+v)*g(l)
#     rdict['history']=[]
#     rdict['logbook']=[]
#     reference_rotation_inner(rdict,p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,0)
# 
# def reference_rotation_inner(rdict,p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,d):
#     maxpow = 20
#     if m >= maxpow:
#         return
# 
#     rdict['history'].append((p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,d))
#     # p,ff,gg,u0,v0,l0,ld,md,m,k,scale,phi,d=rdict['history'][-1]
# 
#     callnumber=len(rdict['history'])-1
#     prefix = " " * d
#     msg=prefix+"Entering recursive call number %d, depth %d" % (callnumber,d)
#     rdict['logbook'].append(msg)
#     print msg
# 
#     a = vector([p^k,0])
#     b = vector([0,p^k])
#     uv0=vector((u0,v0))
# 
#     K=GF(p); KP=K['x']; KR=FractionField(KP)
#     Z=Integers(); ZP=Z['x'];
#     x=ZP.gen()
#     sarr=rdict['sarr']
#     fence=vector([rdict['umax'],rdict['vmax']])
# 
#     h = rdict['h']
#     l,u,v=h.parent().gens()
#     nh = h(phi(l),u0+p^k*u,v0+p^k*v)
#     assert valuation(nh.content(),p)==m
#     assert nh/p^m == ff(l)+(u*phi(l)+v)*gg(l)
#     assert valuation(nh.derivative(l).content(),p)>=m
#     # assert valuation(nh.derivative(l).content(),p)==ld+md
# 
#     # Our expression is c0(l)+u*cu(l)+v*cv(l)
#     c0 = ZP([c % p^maxpow for c in ff.coeffs()])
#     cv = ZP([c % p^maxpow for c in gg.coeffs()])
#     cu = phi * cv
#     assert valuation((nh - p^m*(c0(l)+u*cu(l)+v*cv(l))).content(),p) >= maxpow
#     assert cv == Integers(p^maxpow)['x'](rdict['g'](phi))
# 
#     # The assertions on the trivariate polynomial nh above contain more
#     # that this, but here follows a breakdown just in case.
#     assert p^m*ff(l) == rdict['f'](phi(l))+((u0*phi(l)+v0)*rdict['g'](phi(l)))
#     assert p^m*(u*phi(l)+v)*gg(l)==p^k*((u*phi(l)+v)*rdict['g'](phi(l)))
#     assert p^m*gg(l)==p^k*rdict['g'](phi(l))
# 
#     cc0 = valuation(c0.content(),p)
#     ccu = valuation(cu.content(),p)
#     ccv = valuation(cv.content(),p)
#     assert ccu >= ccv
# 
#     c = min(cc0,ccu,ccv)
#     assert c == 0
# 
#     # zview(mround(matrix(sbound,sbound,sarr)-complete),u0,v0,p^k)
# 
#     dc0 = c0.derivative()
#     dcu = cu.derivative()
#     dcv = cv.derivative()
# 
#     # Now we know that the minimum valuation is zero.
# 
#     g0 = KP(c0)
#     gu = KP(cu)
#     gv = KP(cv)
# 
#     if not g0.is_constant():
#         # Then it's the typical case, as encountered for instance at the
#         # beginning of the root sieve.
#         # Each possible l value must be tried, and will give rise to
#         # potentially intersecting lattices. 
#         for l in range(p):
#             if k==0: l0=l
#             g0l = g0(l); dg0l = KP(dc0)(l)
#             gul = gu(l); dgul = KP(dcu)(l)
#             gvl = gv(l); dgvl = KP(dcv)(l)
#             if gvl == 0: continue
#             msg=prefix+"Looking at roots phi(%r), where phi==%r" %(l,phi)
#             rdict['logbook'].append(msg)
#             print msg
#             igl = 1/gvl
#             du=0
#             dv0=Z(-g0l * igl)
#             dv=dv0
#             u1 = u0 + p^k * du
#             v1 = v0 + p^k * dv
#             uv1=vector((u1,v1))
#             rdict['logbook'].append((uv1, a-phi(l)*b, p* b, -scale/(p-1)))
#             hits = light_array(sarr, uv1, a-phi(l)*b, p* b, fence,-scale/(p-1))
#             msg=prefix+ "%d hits" %hits
#             rdict['logbook'].append(msg)
#             print msg
#             if hits == 0: continue
#             # We're optimistic here ; we're counting on the fact that we
#             # expect the derivative not to cancel.
#             
#             # However, there is the possibility that we reach the
#             # cancellation point. This only once per (p,p) square,
#             # meaning that u,v are constrained.
# 
#             msg=prefix+ "Looking for multiple roots at phi(%r), where phi==%r" %(l,phi)
#             rdict['logbook'].append(msg)
#             print msg
#             # The good equation for u is given by:
#             df=rdict['f'].derivative()
#             dg=rdict['g'].derivative()
# 
#             # Try to guess the right value for du.
#             Zm1=Integers(p^(m+1))
#             # rhs = Z(Zm1((f*dg-df*g)(phi(l))/cv(l)^2-u0))
#             rhs = Z(K(((f*dg-g*df)(phi(l))/cv(l)^2-u0)/p^(m-ld)))
#             # va = valuation(rhs,p)
#             # assert va >= md
#             # assert k >= md
#             # if va < k:
#                 # # Then by no means the derivative will vanish.
#                 # continue
#             range_du=[]
#             if k-m+ld > 0 and rhs == 0:
#                 # Then all u's will do.
#                 range_du=range(p)
#             else:
#                 range_du = [rhs]
#             ## u = Z(-igl*K((((df+(u0*x+v0)*dg+u0*g)(phi))(l)-ff(l)*(p^k*dg(phi))(l))/p^m))
#             ## du = Z((-dg0l * gvl + g0l * dgvl) * igl ^ 2)
#             ## u = Z(Z(Integers(p^(m+1))((f*dg-df*g)(phi(l))/cv(l)^2)-uv0[0])/p^k)
#             for du in range_du:
#                 dv = (dv0 - du *l0) % p
#                 u1 = u0 + p^k * du
#                 v1 = v0 + p^k * dv
#                 uv1=vector((u1,v1))
#                 rdict['logbook'].append((uv1, p*a, p*b, scale/(p-1)/p))
#                 hits = light_array(sarr,
#                             uv1, p*a, p*b, fence, scale/(p-1)/p)
#                 msg=prefix+ "%d hits" %hits
#                 rdict['logbook'].append(msg)
#                 print msg
#                 if hits == 0:
#                     continue
#                 # scale/(p-1)/p is scale/(p-1)-scale/p ; hence it brings
#                 # back the cells to the contribution -scale/p, which is the
#                 # contribution from a _single_ zero at this location, not
#                 # liftable because of ramification. Later recursive calls
#                 # will investigate the possibility that despite the multiple
#                 # root mod p, we still get roots at higher precision.
#                 nphi=phi(l+p*x)
#                 nff=ZP((ff(l+p*x)+(du*nphi+dv)*gg(l+p*x))/p)
#                 ngg=gg(l+p*x)
#                 ## l,u,v=h.parent().gens()
#                 ## nnh=nh(nphi(l),uv1[0]+p*u,uv1[1]+p*v)
#                 ## assert nnh == p*(nff(l)+(u*nphi(l)+v)*ngg(l))
#                 # ff,gg,uv0,l0,m,k,scale,phi=nff,ngg,uv1,l0,m+1,k+1,scale/p,nphi
#                 reference_rotation_inner(rdict,p,nff,ngg,u1,v1,l0,ld+1,md+1,m+1,k+1,scale/p,nphi,d+1)
#     elif ccv > 0:
#         # ccv is zero, so ccu is even larger, which implies that cc0 has
#         # to be zero. Since c0 is a constant polynomial, this implies
#         # that there is no solution. So we finish here.
#         return
#     else:
#         msg=prefix+ "We are here in recursive call number %d" % callnumber
#         rdict['logbook'].append(msg)
#         print msg
#         # cc0 is not zero, which means that our expression has an extra
#         # root only for specific u,v pairs. Basically, we have a linear
#         # term in u and v which must be cancelled.
#         constant_term=g0.constant_coefficient()
#         dv0=Z(-constant_term/gv(l0))
#         dv=dv0
#         u1 = u0
#         v1 = v0 + p^k * dv
#         uv1=vector((u1,v1))
#         assert ccv==0
#         # Furthermore, our expression is c0(l)+u*cu(l)+v*cv(l) ; We know
#         # that cu == phi * cv, hence our equation is simply v == -phi * u
#         msg=prefix+ "Looking for places where the multiple root at %r may lift"%phi
#         rdict['logbook'].append(msg)
#         print msg
#         rdict['logbook'].append((uv1, a-l0*b, p*b,-scale))
#         hits=light_array(sarr, uv1, a-l0*b, p* b, fence,-scale)
#         msg=prefix+ "%d hits" %hits
#         rdict['logbook'].append(msg)
#         print msg
#         if hits == 0:
#             return
#         msg=prefix+"%r"%(K['l','u','v'](c0(l)+u*cu(l)+v*cv(l)))
#         rdict['logbook'].append(msg)
#         print msg
#         for du in range(p):
#             dv = (dv0-du*l0) % p
#             uvscal=vector((du,dv))
#             msg=prefix+ "Looking for lift specifically at u,v === %r" % uvscal
#             rdict['logbook'].append(msg)
#             print msg
#             u2 = u0 + p^k * du
#             v2 = v0 + p^k * dv
#             uv2=vector((u2,v2))
#             msg=prefix+ "global sub-lattice %r+%r,%r" % (uv2,p*a,p*b)
#             rdict['logbook'].append(msg)
#             print msg
#             # We know that phi === l0 mod p, so (u0*phi+v0) is zero mod
#             # p. ff itself is also zero mod p. However, some subtleties
#             # can occur in the difference between the two.
#             nff=ZP((ff+(du*phi+dv)*gg)/p)
#             l,u,v=h.parent().gens()
#             assert nh(l,du+p*u,dv+p*v)/p^(m+1) == nff(l)+(u*phi(l)+v)*gg(l)
#             ## ff,u0,v0,m,k=nff,u1,v1,m+1,k+1
#             reference_rotation_inner(rdict,p,nff,gg,u2,v2,l0,ld,md,m+1,k+1,scale,phi,d+1)
#     msg=prefix+ "Finishing recursive call number %d" % callnumber
#     rdict['logbook'].append(msg)
#     print msg
# 
# 
# def filter_logbook(rdict,u,v):
#     vec=vector((u,v))
#     for s in rdict['logbook']:
#         if type(s) == type(''):
#             print s
#         else:
#             st="sub-lattice %r+%r,%r: %f" % s
#             if ((vec-s[0])*matrix([s[1],s[2]])^-1) in vec.parent():
#                 st="** " + st
#             print st
# 
# def filter_logbook_quick(rdict,u,v):
#     vec=vector((u,v))
#     for s in rdict['logbook']:
#         if type(s) != type(''):
#             if ((vec-s[0])*matrix([s[1],s[2]])^-1) in vec.parent():
#                 st="sub-lattice %r+%r,%r: %f" % s
#                 st="** " + st
#                 print st
# 
# 
# def light_array(sarr,x0,a,b,fence,value):
#     x = copy(x0)
#     pos = x[0] * fence[1]
#     assert b[0] == 0
#     dpos = a[0] * fence[1]
#     hits=0
#     while x[0] < fence[0]:
#         # reduce within the row, since we know that b[0] is zero
#         x[1] = x[1] % b[1]
#         while x[1] < fence[1]:
#             sarr[pos + x[1]] += value
#             hits+=1
#             x[1] += b[1]
#         x += a
#         pos += dpos
#     return hits
# 
# # These are helpers, to be removed once the code works ok.
# def zview(m,i0,j0,d):
#     ncols = ((m.ncols() - j0) / d).ceil()
#     nrows = ((m.nrows() - i0) / d).ceil()
#     n=matrix(nrows,ncols)
#     for i in range(0,nrows):
#         for j in range(0,ncols):
#             n[i,j]=m[i0+i*d,j0+j*d]
#     return n

def mround(m):
    d=matrix(ZZ,m.nrows(),m.ncols())
    for k in range(m.nrows()):
        for l in range(m.nrows()):
            d[k,l]=ZZ(floor(m[k,l]*1000))
    return d

def mprint(m):
    print m.str()

def printmin(arr,X,Y,c):
    min = 1
    u = v = 0
    for x in range(X):
        for y in range(Y):
            if arr[x*Y+y]<min:
                min = arr[x*Y+y]
                u = x
                v = y
    return min+c,u,v

def compose_reduce(f,phi,pmax):
    return f.parent()([c % pmax for c in f(phi).coeffs()])


# Now we work on a trimmed-down version.
def rotation_handle_p(rdict,p):
    """Fast root sieve modulo p and its powers. The rdict argument must
    have been setup with rotation_init beforehand"""
    u0=v0=0
    ff,gg=rdict['f'],rdict['g']
    f,g=rdict['f'],rdict['g']
    rdict['common_contribution']+=float(log(p))/(p-1)
    rdict['common_contribution']+=alpha_p_projective_nodisc(f,p)
    x=Integers()['x'].gen()
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
    rdict['cutoff']=1.0e-5
    maxpow=20
    rdict['maxpow']=maxpow    
    rdict['ppow']=[1,p]
    K=GF(p); KP=K['x'];
    rdict['gmodp']=KP(g)
    fdg_gdf=rdict['fdg_gdf']    
    extend_ppow_list(rdict,maxpow)
    ppow=rdict['ppow']
    fdg_gdf = compose_reduce(fdg_gdf,x,ppow[maxpow])
    ff = compose_reduce(ff,x,ppow[maxpow])
    gg = compose_reduce(gg,x,ppow[maxpow])
    # dphi is (phi-l0)/p -- so this starting value is rubbish, just to
    # get the universe right...
    dphi=x
    rdict['l0']=0
    hits=rotation_inner(rdict,p,ff,gg,fdg_gdf,u0,v0,l0,ld,m,0,scale,dphi,[])
    best,u,v=printmin(rdict['sarr'],rdict['umax'],rdict['vmax'],rdict['common_contribution'])    
    return hits,best,u,v



def rotation_inner(rdict,p,ff,gg,fdg_gdf,u0,v0,l0,ld,m,twist_v1,scale,dphi,hist):
    """Internal routine for the fast root sieve. We are looking at roots
    of the form phi(x), where phi has the form constant_term+p^ld*x. ff,
    gg, and fdg_gdf are basically f composed with phi, taking into
    account the necessary adjustments (see the more extensive code for
    some asserts). dphi and twist_v1 are variables valid only below
    recursion depth 1. m is the recursion depth."""

    #print ("m=%d,ld=%d"%(m,ld))
    #print ("\n--------------------------------------------------------")
    #print ("\n p=%d, m=%d, ld=%d, scale=%2.8f, u0=%d, v0=%d, l0=%d, twist_v1=%d" %(p,m,ld,scale,u0,v0,l0,twist_v1))
    #print ("\nff=%s"%ff)
    #print ("\ngg%s"%gg)
    #print ("\nfdg_gdf=%s"%fdg_gdf)
    #print ("\ndphi=%s"%dphi)
    
    if len(hist) > 6:
        nhist=copy(hist)
        nhist.append([dphi,u0,v0,dphi(0),ld,m])
        # print "Found tree of height %d" % len(nhist)
        # for hi in nhist:
            # print hi


    if scale < rdict['cutoff']:
        return 0

    maxpow=rdict['maxpow']
    ppow=rdict['ppow']
    pm=ppow[m]
    pm1=ppow[m+1]
    pmax=ppow[maxpow]

    #print ("\nmaxpow=%d,pmax=%s"%(maxpow,pmax))

    K=GF(p); KP=K['x'];
    Z=Integers(); ZP=Z['x'];
    x=ZP.gen()
    sarr=rdict['sarr']
    um,vm=rdict['umax'],rdict['vmax']

    # gg is never scaled down, and it is always evaluated at points where
    # the value is not zero mod p. So the valuation of g mod p never
    # vanishes.
    
    # This is always useful.
    minus_f_over_g_modp = KP(ff)
    
    scale1 = scale/(p-1)
    scale2 = scale/p
    scale3 = scale1-scale2

    thits = 0

    if m > 0:
        l0=rdict['l0']
        igl=rdict['igl']
        minus_f_over_g_modp*=-igl
        #print ("\nigl=%d"%igl)
    
    #print ("\nminus_f_over_g_modp=%s"%minus_f_over_g_modp)

    # Deciding whether we're in the general case or not has some
    # subtleties, unfortunately.
    look_many_roots = not minus_f_over_g_modp.is_constant() or   \
                    m == 0 and not rdict['gmodp'].is_constant()

            
    #print ("\nlook_many_roots=%d"%(look_many_roots))

    nhist = copy(hist)
    if not look_many_roots: nhist.append([dphi,u0,v0,l0+p*dphi(0),ld,m])

    if look_many_roots:
        # Then it's the typical case, as encountered for instance at the
        # beginning of the root sieve.
        # Each possible l value must be tried, and will give rise to
        # potentially intersecting lattices. 
        for l in range(p):
            nhist=copy(hist)
            nhist.append([dphi,u0,v0,l0+p*dphi(l),ld,m])
            minus_f_over_g_modp_l = minus_f_over_g_modp(l);
            if m==0:
                l0=l
                gvl = rdict['gmodp'](l);
                if gvl == 0: continue
                igl = 1/gvl
                rdict['l0']=l0
                rdict['igl']=igl
                twist_v1 = - l0
                minus_f_over_g_modp_l *= -igl

            #print("\nminus_f_over_g_modp_l=%d"%minus_f_over_g_modp_l)

            dv0=Z(minus_f_over_g_modp_l)
            u1 = u0
            v1 = v0 + pm * dv0

            hits = light(sarr, u1,v1, pm, pm1, twist_v1, um,vm,-scale1)

            # print ("1:m=%d ld=%d hits>0=%d u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f"%(m,ld,(1 if (hits>0) else 0),u1,v1,pm,pm1,twist_v1,-scale1))
            #print("1:%d "%(1 if (hits>0) else 0))


            if hits == 0: continue
            # print "main, level %d : %d hits" % (m, hits)
            thits += hits
            # We're optimistic here ; we're counting on the fact that we
            # expect the derivative not to cancel.
            
            # However, there is the possibility that we reach the
            # cancellation point. This only once per (p,p) square,
            # meaning that u,v are constrained.

            # Try to guess the right value for du.
            # Rp = Integers(ppow[m-ld+1])
            # RpP=Rp['w']
            # rhs = Z(Z((RpP(fdg_gdf) - u0*RpP(gg)^2)(l))/ppow[m-ld])
            # Slightly faster when not doing reductions while evaluating.
            # Most probably due to the cost of computing remainders all
            # the way to the final value, which ends up being more
            # expensive than it should.

            # TODO: There has to be a way to remove valuations very early
            # on here. In fact, my guess is that it could very probably
            # be done every time m-ld increases...
            rhs = Z(K((fdg_gdf(l)-u0*gg(l)^2)/ppow[m-ld]))
            #print ("\nrhs=%d"%rhs)
            if ld > 0 and rhs != 0:                
                continue           

            psi=l+p*x
            nfdg_gdf = compose_reduce(fdg_gdf,psi,pmax)
            
            #print("\nnfdg_gdf=%s"%nfdg_gdf)
            
            nff = compose_reduce(ff,psi,pmax)
            
            #print("\nnff=%s"%nff)
            
            ngg = compose_reduce(gg,psi,pmax)
            
            #print("\nngg=%s"%ngg)
            

            if m == 0:
                # XXX This is fairly bizarre. Results are correct with
                # this, but I don't see at the moment why we don't have
                # psi instead of x here.
                ndphi = x
            else:
                ndphi = dphi(psi)

            
            ##print("ndphi is \n %s"%ndphi)
            

            nff=ZP((nff+dv0*ngg)/p)
            
            ##print("nff in the lattice is %s"%nff)
            

            twist_f=ndphi*ngg
            
            ##print("twisted_f is %s"%twist_f)
            

            if ld > 0:
                if rhs == 0:
                    # Then all u's will do. Note that this case happens
                    # only rather deep inside the recursion. So indeed
                    # we're going to re-hit the places which have been
                    # hit before, but that's not too worrying.
                    nspots=p
            else:
                du0 = Z(rhs*igl^2)
                u1 += du0 * pm
                v1 += du0 * twist_v1
                nff+= du0 * twist_f
                
                ##print("du0=%d, u1=%d, v1=%d"%(du0,u1,v1))
                ##print("nff after twist ld=0 is \n%s"%nff)
                
                nspots=1

            new_twist_v1 = p * twist_v1

            for du in range(nspots):
                hits = light(sarr, u1,v1, pm1, pm1, 0, um,vm, scale3)
                # print ("2:m=%d ld=%d hits>0=%d u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f"%(m,ld,(1 if (hits>0) else 0),u1,v1,pm1,pm1,0,scale3))
                #print("2:%d "%(1 if (hits>0) else 0))
                # #print "main/excep, level %d : %d hits" % (m, hits)
                thits += hits
                # scale3 is scale/(p-1)-scale/p ; hence it brings
                # back the cells to the contribution -scale/p, which is the
                # contribution from a _single_ zero at this location, not
                # liftable because of ramification. Later recursive calls
                # will investigate the possibility that despite the multiple
                # root mod p, we still get roots at higher precision.
                if hits >0:
                    thits += rotation_inner(rdict,p,nff,ngg,nfdg_gdf,u1,v1,l0,ld+1,m+1,new_twist_v1,scale2,ndphi,nhist)
                nff+=twist_f
                
                ##print("nff after twist du=%d is \n%s"%(du,nff))
                
                u1 += pm
                v1 += twist_v1
    else:
        # f mod p is a constant polynomial.  This means that our
        # expression has an extra root only for specific u,v pairs.
        # Basically, we have a linear term in u and v which must be
        # cancelled.
        assert m>0
        #print("\n-f/g mod p is constant")
        dv0=Z(minus_f_over_g_modp.constant_coefficient())
        u1 = u0
        v1 = v0 + pm * dv0
        hits=light(sarr, u1, v1, pm, pm1, twist_v1, um, vm, -scale)
        # print ("3:m=%d ld=%d hits>0=%d u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f"%(m,ld,(1 if (hits>0) else 0),u1,v1,pm,pm1,twist_v1,-scale))
        #print ("3:%d "%(1 if (hits>0) else 0))
        if hits == 0:
            return thits
        # #print "secondary, level %d : %d hits" % (m, hits)
        thits += hits
        nff=ZP((ff+dv0*gg)/p)
        
        ##print("nff after primitivization is \n%s"%nff)
        
        twist_f=dphi*gg
        new_twist_v1 = p * twist_v1
        for du in range(p):
            hits = rotation_inner(rdict,p,nff,gg,fdg_gdf,u1,v1,l0,ld,m+1,new_twist_v1,scale,dphi,nhist)
            # #print "secondary/excep, level %d : %d hits" % (m, hits)
            thits += hits
            u1 += pm
            v1 += twist_v1
            nff += twist_f
            
            ##print("nff after twist du=%d is \n%s"%(du,nff))
        
    return thits

def light(sarr,u0,v0,us,vs,skew,um,vm,value):
    """Internal routine for the fast root sieve."""
    u = u0
    v = v0
    pos = u * vm
    dpos = us * vm
    hits=0
    oversize = vm + skew
    cq,cr=oversize.quo_rem(vs)
    carriage_return=oversize-cr
    v = v % vs
    while u < um:
        # reduce within the row, since we know that b[0] is zero
        # v = v % vs
        while v < vm:
            sarr[pos + v] += value         
            hits+=1
            v += vs
        u += us
        v += skew
        pos += dpos
        v -= carriage_return
        assert cr <= v and v < cr + vs
        if v >= vs:
            v -= vs
    return hits


