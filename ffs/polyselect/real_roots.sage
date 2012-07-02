# Usage: count_real_roots(f,true) for number of real roots.
# Option true stays for "with multiplicities", i.e. a double root counts 2.
# The option count_real_roots(f,false) has both a different algorithm and 
# output: it is a slower and counts the total number of polynomials r such
# that for a/b \approx r the degree of N(a,b) is less than its raw value,
# i.e. deg_t(f)+max(deg(b),deg(a))^deg_x(f).




def hensel_lift(f,r,prec,R=GF(2)['t,x'],df=None):
    t,x=R.gens()
    F=R.base_ring()
    At.<t>=F['t']
    S.<x>=At['x']
    if df==None:
        df=S(f).derivative()
    assert f(r) % t^prec == 0
    assert df(r) % t != 0
    r_new=r-f(r)*At(df(r)).inverse_mod(At(t)^(2*prec))
    assert (r-r_new) % t^prec == 0
    assert f(r_new) % t^(2*prec) ==0
    return r_new

def roots_mod_t(f):
    Fq=f.base_ring().base_ring()
    R.<t,x>=Fq['t,x']
    U.<x>=Fq['x']
    g=U(R(f)(0,x))
    if g==0:
        return Fq.list()
    return [e[0] for e in g.roots()]


def naive_lift(f,r,prec,R):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    assert f(r) % t^prec == 0
    #fcont=gcd(f.coefficients())
    #v=valuation(fcont,t)
    v=prec
    S=f.parent()
    f_new=S(f(r+ t^v*x)/t^v)
    # using prec instead of v loses (v-prec) iterations of naive_lift
    final=[]
    for e in roots_mod_t(f_new):
        r_new=r+t^v * e
        assert f(r_new) % t^(v+1) == 0
        final.append(r_new)
    return final

def all_roots(f,prec,R):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    df=f.derivative()
    l=roots_mod_t(S(f))
    if prec == 1:
        return l
    final=[]
    for r in l:
        if A(df(r)) % t != 0:
            pr=1
            while pr<prec:
                r=hensel_lift(f,r,pr)
                pr*=2
            assert f(r) % t^prec == 0
            final.append(r.truncate(prec))
        else:
            if A(df(r)) % t^2 ==0:
                f_aux=f(t,r+t*x)
                final+=[r+t*e for e in all_roots(f_aux,prec-1,R)]
                continue
            pr=1
            lr=[r]
            while pr<prec:
                lr_new=[]
                for rr in lr:
                    lr_new+=naive_lift(f,rr,pr,R)
                pr+=1
                lr=lr_new
            final+=lr
    return final


def all_roots_rec(f,prec,R=GF(2)['t,x'],pr=None):
    if pr == None:
        pr=prec
    if pr <=0:
        return []
    F=R.base_ring()
    t,x=R.gens()
    At.<t>=F['t']
    S.<x>=At['x']
    f=S(f)
    df=f.derivative()
    l=roots_mod_t(S(f))
    if pr == 1:
        return l
    final=[]
    for r in l:
        if At(df(r)) % t != 0:
            pr0=1
            while pr0<prec:
                r=hensel_lift(f,r,pr0)
                pr0*=2
            assert f(r) % t^prec == 0
            final.append(r.truncate(prec))
        else:
            f_aux=S(f)(r+t*x)
            v_content=valuation(gcd(f_aux.coefficients()),At(t))
            final+=[r+t*e for e in all_roots_rec(f_aux,prec-v_content,R,pr-1)]
    return final

def count_all_roots(f,prec,R):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    S.<x>=A['x']
    f=S(f)
    df=f.derivative()
    l=roots_mod_t(S(f))
    result=0
    for r in l:
        if A(df(A(r))) % A(t) != 0:
            result+=1
        else:
            pr=1
            lr=[r]
            while pr<prec:
                lr_new=[]
                for rr in lr:
                    if gcd(A(df(A(rr))),A(t^pr)) % A(t) != 0:
                        result+=1
                    else:
                        lr_new+=naive_lift(f,rr,pr,R)
                pr+=1
                lr=lr_new
            result+=len(lr)
            # We assume that all the roots at precision prec can be
            # expanded as power series roots
    return result
# given g in Fq(t)[x], rotate(g) is in Fq[t][x] and has its leading
# coefficient of null valuation;
# if r is a root of rotate(g) r*t^s is a root of g

def rotate(g,R):
    F=R.base_ring()
    t,x=R.gens()
    At.<t>=F['t']
    K=FractionField(At)
    S.<x>=K['x']
    g=S(g.numerator())*K(1/g.denominator())
    n,gg=S(g).degree(),S(g).coeffs()
    minv=min([(valuation(gg[i],At(t))-valuation(gg[n],At(t)))/(n-i) for i in range(n)])
    if minv != Infinity:
        s=floor(minv)
    else:
        print "invalid polynomial to rotate"
        assert False
    m=-n*s-valuation(gg[n],At(t))
    g0=t^m*g(x*t^s)
    assert valuation(At(g0.leading_coefficient()),At(t)) == 0
    for c in g0.coefficients():
        assert valuation(At(c),At(t)) >= 0
    g0=R(g0)
    return g0, s, m


def real_roots(f,prec=0,R=GF(2)['t,x']):
    if prec ==0:
        prec=4*f.degree()
    t,x=R.gens()
    F=R.base_ring()
    At.<t>=F['t']
    K=FractionField(At)
    S.<x>=K['x']
    f_bar=R(f)(1/t,x)
    f_tilde,s,m=rotate(f_bar,R)
    roots_f_tilde=find_roots_of_integer_polynomial(f_tilde,R,prec+abs(s),-10000)
    roots_f_bar=[r*t^s for r in roots_f_tilde]
    roots_f=[r(1/t) for r in roots_f_bar]
    return roots_f 


def reconstruct_real_root(r,prec,R):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    K=FractionField(A)
    S.<x>=A['x']
    s=A(K(r).numerator()).degree()-A(K(r).denominator()).degree()
    F1=A(K(r/t^s)(1/t))
    n=ZZ(prec/2)
    M=A(t^prec)
    h=half_xgcd(n,F1,M,A,t)
    # a/b=F1 mod M
    a=h[0]
    b=h[1]
    b*=t^s
    m=max(a.degree(),b.degree())
    a0=A(t^m*a(1/t))
    b0=A(t^m*b(1/t))
    # (a0/b0)(t)=(a/b)(1/t)
    return a0,b0

def test_decrease_in_degree(f,r,prec,a0,b0,R):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    K=FractionField(A)
    S.<x>=A['x']
    df=S(f).derivative()
    d3f=S(f).derivative().derivative().derivative()
    d1=-(df(r).numerator().degree()-df(r).denominator().degree())
    d3=-(d3f(r).numerator().degree()-d3f(r).denominator().degree())
    if 2*prec<=d3-d1:
        print "precision not large enough"
        return False
    if A(b0^f.degree()*f(a0/b0)).degree()>f.degree(x)*b.degree()-prec+d1:
        print "error"
        return False
    return True

#i is an index in range(n+1)
# fc is the list of coefficients in Fq[t] of f
def next_newton_segment(i,fc,A):
    t=A.gen()
    if fc[i] == 0:
        return i-1, -Infinity, 1
    vertex=i-1
    sl=( valuation(fc[i-1],t)- valuation(fc[i],t) )/(i-(i-1))
    for j in range(i):
        sl_new= ( valuation(fc[j],t)- valuation(fc[i],t) )/(i-j)
        if sl_new < sl:
            sl=sl_new
            vertex=j
    return vertex,sl,i-vertex

def count_roots_of_integer_polynomial(f,R,prec,minslope):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    K=FractionField(A)
    S.<x>=A['x']
    if prec==0:
        prec=10*S(f).degree()
    if R(f).degrees()[0] == 0:
        X.<x>=F['x']
        g=X(f)
        return len(g.roots())
    result=0
    fc=S(f).coeffs()
    i=S(f).degree()
    while i>0:
        i,sl,d=next_newton_segment(i,fc,A)
        sl=min(sl,prec)
        if sl in ZZ and sl>= minslope:
            f_new=R(f)(t,x-t^sl)
            if sl < prec:
                result+=count_roots_of_integer_polynomial(f_new,R,prec,sl+1)
            else:
                result+=d
    return result 

def find_roots_of_integer_polynomial(f,R,prec,minslope):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    K=FractionField(A)
    S.<x>=A['x']
    if prec==0:
        prec=10*S(f).degree()
    if R(f).degrees()[0] == 0:
        X.<x>=F['x']
        g=X(f)
        return [e[0] for e in g.roots()]
    result=[]
    fc=S(f).coeffs()
    i=S(f).degree()
    while i>0:
        i,sl,d=next_newton_segment(i,fc,A)
        sl=min(sl,prec)
        if sl in ZZ and sl>= minslope:
            f_new=R(f)(t,x-t^sl)
            if sl < prec:
                result+=[ t^sl+e for e in
                        find_roots_of_integer_polynomial(f_new,R,prec,sl+1)]
            else:
                result+=[0]
    return result 
# If double roots count as 2 roots than multiplicities=True
# and in this case we might speed up computations.
# If two roots are equal by multiplicity t^70 we count as a double root.
# IMPORTANT REMARK: A "root" at precision 10 may not expand to a real root.
# Thus multiplicities=False may find many more roots that the True option
def count_real_roots(f,multiplicities=True,prec=0,R=GF(2)['t,x'],minslope=-100000):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    K=FractionField(A)
    S.<x>=A['x']
    if prec==0:
        prec=10*S(f).degree()
    f_bar=R(f)(1/t,x)
    f_tilde,s,m=rotate(f_bar,R)
    if R(f_tilde).degrees()[0] == 0:
        X.<x>=F['x']
        g=X(f_tilde)
        return len(g.roots())
    if multiplicities == False:
        return count_all_roots(f_tilde,prec,R)
    else:
        return count_roots_of_integer_polynomial(f_tilde,R,prec,-100000)      


def construct_f_with_real_roots(P=0,R=GF(2)['t,x']):
    F=R.base_ring()
    t,x=R.gens()
    A.<t>=F['t']
    K=FractionField(A)
    S.<x>=K['x']
    if P == 0:
        P=A.random_element(20)
    f0=1
    for i in F.list():
        f0*=(x-i)
    f=f0+t*P
    f_bar=R(f)(1/t,x)
    f_tilde,s,m=rotate(S(f_bar),R)
    return f_tilde


