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
    #rotation_inner(f,g,range(-U,U+1),range(-V,V+1))


def rotation_init(f,g,U0,U1,V0,V1):
    """Sets up rotation for using the area [u0..u1[ x [v0..v1["""
    assert(U0==0)
    assert(V0==0)
    space=(U1-U0)*(V1-V0)
    #print "Rotation space: %r" % space 
    rdict=dict(umax=U1,vmax=V1,f=f,g=g)
    rdict['U0'] = U0
    rdict['U1'] = U1
    rdict['V0'] = V0
    rdict['V1'] = V1
    # The list of lattices.
    rdict['lattice_list']=[]
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



def mround(m):
    d=matrix(ZZ,m.nrows(),m.ncols())
    for k in range(m.nrows()):
        for l in range(m.nrows()):
            d[k,l]=ZZ(floor(m[k,l]*1000))
    return d

def mprint(m):
    print m.str()

def printmin(arr,X,Y):
    min = 1
    u = v = 0
    for x in range(X):
        for y in range(Y):
            if arr[x*Y+y]<min:
                min = arr[x*Y+y]
                u = x
                v = y
    print "Minimum is %f, pair is %d,%d" %(min,u,v)

def compose_reduce(f,phi,pmax):
    return f.parent()([c % pmax for c in f(phi).coeffs()])


# Now we work on a trimmed-down version.
def rotation_handle_p(rdict,p):
    """Fast root sieve modulo p and its powers. The rdict argument must
    have been setup with rotation_init beforehand"""
    u0=v0=0
    ff,gg=rdict['f'],rdict['g']
    f,g=rdict['f'],rdict['g']
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
    sarr2 = [float(0.0) for i in range(rdict['umax']*rdict['vmax'])]
    #for lattice in rdict['lattice_list']:
    #   light(sarr2,lattice['u0'],lattice['v0'],lattice['us'],lattice['vs'],lattice['skew'],rdict['umax'],rdict['vmax'],lattice['contrib'])
    for lattice in rdict['lattice_list']:
        light_rectangle(sarr2,lattice['u0'],lattice['v0'],lattice['us'],lattice['vs'],lattice['skew'],0,rdict['umax']-1,0,rdict['vmax']-1,lattice['contrib'])
    if len(rdict['sarr'])==len(sarr2):
        print "Arrays are of same length"
    for i in range(rdict['umax']*rdict['vmax']):
        if rdict['sarr'][i]<>sarr2[i]:
            print "Warning : different values of alpha in entry %d"%i
    print "Finished reviewing average exponent values."
    printmin(rdict['sarr'],rdict['umax'],rdict['vmax'])
    #print_lattices(rdict)
    return hits



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
        ## print "Found tree of height %d" % len(nhist)
        ## for hi in nhist:
        ## print hi


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
    U0,U1,V0,V1=rdict['U0'],rdict['U1'],rdict['V0'],rdict['V1']
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
            

            hits = light_rectangle(sarr, u1,v1, pm, pm1, twist_v1, U0,U1-1,V0,V1-1,-scale1)
            # print ("1:m=%d ld=%d hits>0=%d u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f"%(m,ld,(1 if (hits>0) else 0),u1,v1,pm,pm1,twist_v1,-scale1))
            assert ((1 if (hits>0) else 0) == add_lattice(rdict,u1,v1,pm,pm1,twist_v1,0,um-1,0,vm-1,-scale1,p))
            #hits = add_lattice(rdict,u1,v1,pm,pm1,twist_v1,0,um,0,vm,-scale1,p)

            
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
                hits = light_rectangle(sarr, u1,v1, pm1, pm1, 0, U0,U1-1,V0,V1-1, scale3)
                # print ("2:m=%d ld=%d hits>0=%d u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f"%(m,ld,(1 if (hits>0) else 0),u1,v1,pm1,pm1,0,scale3))
                assert (1 if (hits>0) else 0) == add_lattice(rdict,u1,v1,pm1,pm1,0,0,um-1,0,vm-1,scale3,p)                
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
        hits = light_rectangle(sarr, u1, v1, pm, pm1, twist_v1, U0,U1-1,V0,V1-1, -scale)        
        # print ("3:m=%d ld=%d hits>0=%d u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f"%(m,ld,(1 if (hits>0) else 0),u1,v1,pm,pm1,twist_v1,-scale))
        assert (1 if (hits>0) else 0) == add_lattice(rdict,u1,v1,pm,pm1,twist_v1,0,um-1,0,vm-1,-scale,p)        
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


def add_lattice(rdict,u0,v0,us,vs,skew,U0,U1,V0,V1,contrib,p):
    """Internal routine for testing the fast root sieve."""
  
    if(vs<0):
        vs = -vs

    if(us<0):
        us = -us
        skew = -skew
  
    while(skew<0):
        skew += vs

    while(skew>=vs):
        skew-=vs

    assert vs>0 and skew >= 0 and us>0  

    q = floor((u0-U0)/us)
    u = u0-q*us
    v = v0-q*skew

    q = floor((v-V0)/vs);
    v = v - q*vs;
    
    assert v>=V0
    assert u>=U0

    if not empty(u,v,us,vs,skew,U0,U1,V0,V1):
        lattice = dict(p=p,u0=u0,v0=v0,us=us,vs=vs,skew=skew,contrib=contrib)
        rdict['lattice_list'].append(lattice)        
        return 1
    else:
        return 0

def empty(u0,v0,us,vs,skew,U0,U1,V0,V1): #We suppose u>=U0 and v>=V0
    """ Returns 0 if the given rectangle does not contain a point of the lattice, 1 otherwise"""
    if u0<=U1:
        if v0<=V1:
            return 0
        else:
            if skew==0:
                return 1            
            #q1=ceil((V0-(v0-vs))/skew)
            q2=floor((U1-u0)/us)
            u = u0
            v = v0-vs
            for i in range(q2):                
                    u+=us
                    v+=skew
                    q = floor((v-V0)/vs)
                    v = v-q*vs
                    if v<=V1 and u<=U1:
                        #sys.stderr.write("u0=%d v0=%d u=%d v=%d i=%d\n"%(u0,v0,u,v,i))
                        return 0            
            return 1
    else:
        return 1

def print_lattices(rdict):
    for lattice in rdict['lattice_list']:
        print "u0=%d v0=%d us=%d vs=%d skew=%d contrib=%f p=%d"%(lattice['u0'],lattice['v0'],lattice['us'],lattice['vs'],lattice['skew'],lattice['contrib'],lattice['p'])



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

def light_rectangle(sarr,u0,v0,us,vs,skew,U0,U1,V0,V1,value):
    """ Same as light, but in a rectangle. """
    if(vs<0):
        vs = -vs

    if(us<0):
        us = -us
        skew = -skew
  
    while(skew<0):
        skew += vs

    while(skew>=vs):
        skew-=vs

    assert(us>0 and vs>0 and skew>=0)

    q = floor((u0-U0)/us)
    u = u0-q*us
    v = v0-q*skew

    q = floor((v-V0)/vs);
    v = v - q*vs;
    
    assert v>=V0
    assert u>=U0
    
    hits=0    
    pos=(V1-V0+1)*(u-U0) + (v-V0)

    while u<=U1:
        while v<=V1:
            sarr[pos]+=value
            hits+=1
            v+=vs
            pos+=vs
        u+=us
        v+=skew
        q=floor((v-V0)/vs)
        v=v-q*vs
        pos=(V1-V0+1)*(u-U0) + (v-V0)

    return hits


        
