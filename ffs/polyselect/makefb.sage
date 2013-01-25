load l_roots.sage
load hexa.sage

"""
 input : lr a pair of lists containing pairs [r,k] s. t. F(r,1) % l^k=0 and
 F(1,r) % l^k =0 respectively.
 output: 
 simple is a list of roots r such that f(r) % l =0 and f'(r) % l !=0 
 and possibly r=l if F(1,l) % l =0 and F(1,l) % l^2 !=0

 ram is a list in makefb format i.e. each line is 
[l,r]    if F(r,1) % l = 0 and  f'(r) % l !=0
[l, l+r] if F(1,r) % l = 0 and d(F(1,r))/dr !=0
[l^k,r,v1,v0] if val(F(r,1),l)=v1  and  val(F( r % l^(k-1),1),l)=v0 
[l^k,l^k+r,v1,v0] if val(F(1,r),l)=v1  and  val(F(1,r % l^(k-1)),l)=v0 
"""
def makefb_format(f,l,lr):
    A=l.parent();t=A.gen();F=A.base_ring();R.<t,x>=F['t,x'];S.<x>=A['x'];
    f=S(f)
    df=f.derivative()
    cf=f.coeffs()
    cf.reverse()
    f_proj=S(cf) # f_proj(x)=F(t,1,x), f=F(t,x,1)
    df_proj=f_proj.derivative()
    aff=lr[0]
    proj=lr[1]
    ram=[]
    simple=[]
    for r_ in aff:
        r=r_[0]
        k=r_[1]
        if k == 1 and A(df(r)) % A(l) != 0:
            simple.append(r)
        else:
            if (A(r) % A(l)) in simple:
                continue
            v1=valuation(A(S(f)(r)),l)
            if k == 1:
                v0=0
            else:
                v0=valuation(A(S(f)(A(r) % A(l)^(k-1))),l)
            if v1 > v0:
                ram.append([l^k,v1,v0,r]) 
    for r_ in proj:
        r=r_[0]
        k=r_[1]
        if k == 1 and A(df_proj(r)) % A(l) != 0:
            assert r == 0
            simple.append(r+l)
        else:
            if ((A(r) % A(l))+A(l)) in simple:
                continue
            v1=valuation(A(S(f_proj)(r)),l)
            if k == 1:
                v0=0
            else:
                v0=valuation(A(S(f)(r % l^(k-1))),l)
            ram.append([l^k,v1,v0,r+l^k]) 
    return ram,simple

"""
output: a list written one element per line in file filename or on screen
dlim and powerlim are such that items [l^k,...] verify deg(l)<=dlim and
deg(l^k)<=powerlim. It is required that powerlim <= dlim.
"""
def makefb(f,dlim,powerlim,filename=""): 
    A=f.base_ring();t=A.gen();F=A.base_ring();R.<t,x>=F['t,x'];S.<x>=A['x'];
    if filename == "":
        gd=sys.stdout
    else:
        gd=open(filename,"w")
    gd.flush()
    gd.write("# f="+str(f)+"\n")
    gd.write("# dlim="+str(dlim)+"\n")
    gd.write("# powerlim="+str(powerlim)+"\n")
    gd.flush()
    assert powerlim <= dlim
    for l in A.polynomials(max_degree=dlim):
        if not A(l).is_irreducible():
            continue
        prec=max(floor(powerlim/l.degree()),1)
        ram,simple=makefb_format(S(f),A(l),all_roots(S(f),A(l),prec))
        c=len(simple)
        if c > 0:
            output=hexify(l)+": "+hexify(A(simple[0]))
            for i in range(1,c):
                output+=","+hexify(A(simple[i]))
            gd.write(output+"\n")
        for item in ram:
            gd.write(hexify(A(item[0]))+":"+str(item[1])+","+str(item[2])+": "+hexify(A(item[3]))+"\n")
        gd.flush()
    if gd != sys.stdout:
        gd.close()

""" EXAMPLE
makefb_format(S((x+t)^3+t^5),A(t),all_roots(S(x+t),A(t),7))
makefb(S(t*x^3-(1+t)^2),5,3)
"""
