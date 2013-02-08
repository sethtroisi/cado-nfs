# tools.sage contains functions related to representation of polynomials
# (e.g hexadecimal), listing of irreducible monic polynomials, reduction of
# a polynomial in 2 variables modulo a polynomial in one variable etc.


# Primes outputs a the list of all the irreducible monic polynomials in A
# such that the representation of each element as an integer is in the
# interval [astart,alim]. 
def Primes(A,alim,astart=0):
	if A==ZZ:
		return primes(astart,alim)
	else:
		res=[]
		for a in element_list(A,alim,astart):
			if (not a.is_zero()) and (a.is_irreducible()) and (a.is_monic()):
				res.append(a)
	return res
    
# bar computes the reduction of f \in A[x] modulo p \in A. 
def bar(f,l):
    A=f.base_ring()
    if A==ZZ:
        res=f.change_ring(GF(l))
    else:
        t=A.gen()
        F=A.base_ring()
        p=F.characteristic()
        kl.<w>=GF(F.cardinality()^l.degree(),name='w')
        ww=F.modulus().roots(kl)[0][0]
        if w==GF(p)(1):
            w=l.roots(GF(l))[0][0] # convention of sage is w=1 when the residual field is GF(q)
        newX.<x>=kl['x']
        res=0
        i=0
        for c in f.coeffs():
            res+=c.change_ring(kl)(t=w)*x^i
            i+=1
        res=newX(res)
    return res

# lon(Norm(p)) as described below.
def log_norm(p):
    A=p.parent()
    if A==ZZ:
        res=log(p)
    else:
        F=A.base_ring()
        q=len(F)
        res=p.degree()*log(q)
    return res

# Norm(p) is the cardinal of A/<p>.
def Norm(p):
    A=p.parent()
    if A==ZZ:
        res=abs(p)
    else:
        F=A.base_ring()
	q=len(F)
        res=q^p.degree()
    return res

# element_list outputs the list of all the elements of A[x] such that their
# integer representation belongs to the interval [astart,alim]
def element_list(A,alim,astart=0):
	if A==ZZ:
		return range(astart,alim)
	else:
		q=len(A.base_ring())
		F=A.base_ring()
		t=A.gen()
		l=F.list()
		res=[]
		for i in range(astart,alim):
			iPol=int2pol(q,ZZ(i).digits(q),l,A)			
			res+=[iPol]
		return res			

# 
def int2pol(q,dig,l,A):
# dig =list of digits of an integer in base q
# l =list of elements of GF(q)
# ring where belongs the result
   res=A(0)
   t=A.gen()
   k=0
   for j in dig:
       res+=l[j]*t^k
       k+=1
   return res


def cast(d,A,q): # GF(q) to ZZ for small q
	l=element_list(A,q)
	i=0
	for j in l:
		if j==d:
			return i
		i+=1
		if i>q:
			j='error'
	print 'cast error'
	return 0	
	
def position(a):
	if (a+0).parent()==ZZ:#int+0=integer
		return a
	else:
		A=(a+0).parent()
		q=ZZ(len(A.base_ring()))
		res=0
		l=a.coeffs()
		k=0
		for d in l:
			res+=cast(d,A,q)*q^k
			k+=1 		
		return res



def parse_line(l,A=ZZ):
	l=l.split('\n')[0] # eliminating '\n'	
	words=l.split(':')
	q=A(words[0])
	if len(words)==3:
		n=[ZZ(i) for i in words[1].split(',')]
		r=[A(ri) for ri in words[2].split(',')] 
	else:
		n=[1,0]
		r=[A(ri) for ri in words[1].split(',')]
	return q,n,r 	

def homogenize(f):
    A=f.parent()
    R.<X,Y>=A['X,Y']
    res=R(0)
    cl=f.coeffs()
    d=f.degree()
    i=0
    for i in range(d+1):
        res+=X^i*Y^(d-i)*cl[i]
    return res
    


#def clearall():
#   all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")] 
#   for var in all:
#        del globals()[var]


def fmodp_roots(fcoeffs,l,K,Ky):
    if len(fcoeffs)==2: # degree 1 polynomial
       den = fcoeffs[1] % l
       if den <> 0:
          return [(-fcoeffs[0] * den.xgcd(l)[1]) % l]
       else: # lc(f) % l == 0:
          return []
    A=l.parent()
    t=A.gen()
    p=A.characteristic()
    F=l.base_ring()
    q=F.cardinality()
    wq=F.gen()
    c=K.cardinality()
    wK=K.gen()
    assert c==q^l.degree()
    y=Ky.gen()
    wqK=F.modulus().roots(K)[0][0]
    lcoeffs=A(l).coeffs()
    lK=Ky(0)
    for j in range(l.degree()+1):
        lj=lcoeffs[j]
        if q == p^1:
            ljK=lj
        else:
            ljK=K(0)
            lj_vector=F(lj)._vector_()
            for h in range(len(lj_vector)):
                ljK+=lj_vector[h]*wqK^h
        lK+=ljK*y^j
    w=lK.roots(K)[0][0]
    f_mod_l=Ky(0)
    for i in range(len(fcoeffs)):
        fi=fcoeffs[i]
        ficoeffs=fi.coeffs()
        if ficoeffs == []:
            ficoeffs=[0]
        fiK=K(0)
        j=0
        fij=ficoeffs[j]
        if q == p^1:
            fijK=fij
        else:
            fij_vector=F(fij)._vector_()
            fijK=K(0)
            for h in range(len(fij_vector)):
                fijK+=fij_vector[h]*wqK^h
        fiK+=fijK*1 # w^j=1 including the case w=0
        for j in range(1,fi.degree()+1):
            fij=ficoeffs[j]
            if q == p^1:
                fijK=fij
            else:
                fij_vector=F(fij)._vector_()
                fijK=K(0)
                for h in range(len(fij_vector)):
                    fijK+=fij_vector[h]*wqK^h
            fiK+=fijK*w^j
        f_mod_l+=fiK*y^i
    roots_mod_l=f_mod_l.roots()
    result=[]
    base=[]
    i=0 # w^i=1 including the case w=0
    for j in range(F.degree()):
        if K.degree() == 1:
            base.append([1*wqK^j])
        else:
            base.append(K(1*wqK^j)._vector_())
    for i in range(1,l.degree()):
        base.append(K(w^i)._vector_())
        for j in range(1,F.degree()):
            if K.degree() == 1:
                base.append([w^i*wqK^j])
            else:
                base.append(K(w^i*wqK^j)._vector_())
    matrix=(Matrix(base).transpose())^(-1)
    for rr in roots_mod_l:
        r=rr[0]
        if K.degree() == 1:
            r_vect=vector([r])
        else:
            r_vect=K(r)._vector_()
        r_wrt_w_vect= matrix*r_vect
        r_wrt_w=sum([sum([GF(p)(r_wrt_w_vect[i*F.degree()+j])*F(wq)^j*A(t)^i for j in range(F.degree())]) for i in range(l.degree())])
        #assert A(S(f)(r_wrt_w)) % A(l) == 0
        result.append(r_wrt_w)
    return result





    res=Ky(0)
    for i in range(len(fcoeffs)):
        cc=A(fcoeffs[i])
        cc_mod_l=K(0)
        cc_coeffs=cc.coeffs()
        cc_0=cc_coeffs[0].polynomial()(wq=wqK)
        cc_mod_l+=cc_0
        for j in range(1,len(cc_coeffs)):
            cc_j=cc_coeffs[j].polynomial()(wq=wqK)
            cc_mod_l+=w^j*cc_j
        res+=cc_mod_l*y^i
    r=res.roots()
    if l.degree()==1:
        return [el[0] for el in r]
    matrix=Matrix([K(w^i)._vector_() for i in range(K.degree())]).transpose()
    ll=list(matrix^(-1)*vector([0,1]+[0]*(K.degree()-2)))
    phi=A(0)   # w0=phi(w)
    for i in range(len(ll)):
        phi+=t^i*ll[i]
    #assert w0==sum(w^i*K(phi.coeffs[i]) for i in range(phi.degree()+1)) 
    r_pol_format=[]
    for el in r:
        ri=el[0]
        g=A(K(ri).polynomial())
        r_pol_format.append(g(phi(t)) % l)
    return r_pol_format



def sigma(hexpol,q=2,A=None):
    if A==None:
        A.<t>=GF(q)['t']
    else:
        t=A.gen()
    l=A.base_ring().list()
    dummy=2^ceil(log(q,2))
    return int2pol(dummy,ZZ(int(hexpol,16)).digits(dummy),l,A)


def hexify(f):
    p=f.base_ring().characteristic()
    A.<t>=GF(p)['t']
    S.<x>=A['x']
    fc=S(f).coeffs()
    dummy=2^ceil(log(p)/log(2))
    result=''
    aux=0
    if f.degree() < 1:
        return hex(ZZ(f))
    l=fc[0].coeffs()
    for j in range(len(l)):
        aux+=ZZ(l[j])*dummy^j
    result=result+hex(aux)
    for i in range(1,len(fc)):
        aux=0
        l=fc[i].coeffs()
        for j in range(len(l)):
            aux+=ZZ(l[j])*dummy^j
        result=result+','+hex(aux)
    return result

