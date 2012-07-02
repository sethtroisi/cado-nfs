
def init_rings(p):
    F=GF(p)
    A.<t>=F['t']
    K.<t>=FractionField(A)
    S.<x>=A['x']
    R.<t,x>=F['t,x']
    return F,A,K,S,R


def Primes(A,alim,astart=0):
	if A==ZZ:
		return primes(astart,alim)
	else:
		res=[]
		for a in element_list(A,alim,astart):
			if (not a.is_zero()) and (a.is_irreducible()) and (a.is_monic()):
				res.append(a)
	return res
    

def next_irreducible(P):
    A=P.parent()
    F=A.base_ring()
    q=F.cardinality()
    Pint=position(P) # representation of P as a number
    Qint=Pint+1
    Q=int2pol(q ,ZZ(Qint).digits(q),F.list(),A)
    while not Q.is_irreducible():
        Qint+=1
        Q=int2pol(q,ZZ(Qint).digits(q),F.list(),A)   
    return Q

def digits(phi,m):
    A=phi.parent()
    result=[]
    while phi !=0 :
        quo,rem=phi.quo_rem(m)
        result.append(rem)
        phi= quo
    if result == []:
        result=[0]
    return result

def bar(f,p):
    A=f.base_ring()
    if A==ZZ:
        res=f.change_ring(GF(p))
    else:
        t=A.gen()
        q=A.base_ring().characteristic()
        kp.<w>=GF(q^p.degree(),name='w',modulus=p)
        if w==GF(q)(1):
            w=p.roots(GF(q))[0][0] # convention of sage is w=1 when the residual field is GF(q)
        newX.<x>=kp['x']
        res=0
        i=0
        for c in f.coeffs():
            res+=c.change_ring(kp)(t=w)*x^i
            i+=1
        res=newX(res)
    return res


def log_norm(p):
    A=p.parent()
    if A==ZZ:
        res=log(p)
    else:
        F=A.base_ring()
	q=len(F)
        res=p.degree()*log(q)
    return res

def Norm(p):
    A=p.parent()
    if A==ZZ:
        res=abs(p)
    else:
        F=A.base_ring()
	q=len(F)
        res=q^p.degree()
    return res

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

def parse_line(l,A=ZZ,typo="human"):
    l=l.split('\n')[0] # eliminating '\n'	
    words=l.split(':')
    if typo == 'hex':
        pp=sigma(words[0],p,A)
    else:
        pp=A(words[0])
    if len(words)==3:
        n=[ZZ(i) for i in words[1].split(',')]
        if typo == 'hex':
            r=[sigma(ri,p,A) for ri in words[2].split(',')] 
        else: 
            r=[A(ri) for ri in words[2].split(',')] 
    else:
        n=[1,0]
        if typo == 'hex':
            r=[sigma(ri,p,A) for ri in words[1].split(',')]
        else:
            r=[A(ri) for ri in words[1].split(',')]
	return pp,n,r 	

def homogenize(f):
    A=f.base_ring()
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


def fmodp_roots(fcoeffs,p,K,Ky):
    if len(fcoeffs)==2: # degree 1 polynomial
       den = fcoeffs[1] % p
       if den <> 0:
          return [(-fcoeffs[0] * den.xgcd(p)[1]) % p]
       else: # lc(f) % p == 0:
          return []
    A=p.parent()
    t=A.gen()
    q=p.base_ring().cardinality()
    c=K.cardinality()
    w0=K.gen()
    assert c==q^p.degree()
    y=Ky.gen()
    w=p.roots(K)[0][0]
    l=[c(t=w) for c in fcoeffs]
    res=Ky(0)
    for i in range(len(l)):
        res+=l[i]*y^i
    r=res.roots()
    if p.degree()==1:
        return [el[0] for el in r]
    matrix=Matrix([K(w^i)._vector_() for i in range(p.degree())]).transpose()
    l=list(matrix^(-1)*vector([0,1]+[0]*(p.degree()-2)))
    phi=A(0)   # w0=phi(w)
    for i in range(len(l)):
        phi+=t^i*l[i]
    assert w0==phi.base_extend(K)(w) 
    r_pol_format=[]
    for el in r:
        ri=el[0]
        g=A(K(ri).polynomial())
        r_pol_format.append(g(phi(t)) % p)
    return r_pol_format
# S([sigma(e) for e in hexify(f).split(',')]) == f
def sigma(hexpol,q=2,A=None):
    if A==None:
        A.<t>=GF(q)['t']
    else:
        t=A.gen()
    l=A.base_ring().list()
    return int2pol(2,ZZ(int(hexpol,16)).digits(2),l,A)

def random_multivariate_polynomial(R,degrees_list,terms,sl=0):
    if len(R.gens()) != len(degrees_list):
        print "incorrect number of degrees in random_multivariate_polynomial"
        return 0
    result=R(0)
    for i in range(terms):
        term=R(1)
        for j in range(len(degrees_list)):
            nj=randint(0,degrees_list[j])
            term*=R.gens()[j]^nj
        result+=term
    return result

def random_univariate_polynomial(S,d,dt,sl=0):
    # S = ring of polynomials
    # d = degree 
    # dt = coefficient size
    # sl = slope i.e. result[i].degree(t)+i*sl<dt
    x=S.gen()
    A=S.base_ring()
    result=S(0)
    for i in range(d+1):
        if A == ZZ:
            coeff=randint(0,dt)
        else:
            coeff=A.random_element(max(dt-sl*i,0))
        result+=x^i*coeff
    return result


def hexify(f):
    p=f.base_ring().characteristic()
    A.<t>=GF(p)['t']
    S.<x>=A['x']
    fc=S(f).coeffs()
    dummy=2^ceil(log(p)/log(2))
    result=''
    aux=0
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
