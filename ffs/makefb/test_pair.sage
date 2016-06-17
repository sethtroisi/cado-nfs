load tools.sage
load makefb.sage

def test_pair(ab_list,powerlim,F,filename):
	N=len(ab_list)
	A=ab_list[0][0].parent()
	val=[0 for i in range(N)]
	fd=open(filename,"r").readlines()
	p=1
	unhit=0
	total=0
	for line in fd:
		if line[0]=='#':
			continue 	
		pk,n,r=parse_line(line,A)
		spy=true
		if (p ==1):
			p=pk.factor()[0][0]
			kmax=ceil(powerlim*log(2)/log_norm(p))
			spy=false
		if spy and (pk % p !=0):
			po=p
			kmax=ceil(powerlim*log(2)/log_norm(p))
			p=pk.factor()[0][0]
			for i in range(N):
				a,b=ab_list[i][0],ab_list[i][1]
				if gcd(a,b)==1:
					total+=1
					if F(a,b)!=0 and min(valuation(F(a,b),po),kmax-1) != val[i]:
						print "F,a,b,p",F,a,b,po,val[i],valuation(F(a,b),po),kmax
						unhit+=1
			val=[0 for i in range(N)]
		for ri in r:
			if position(ri)<position(pk):
				ua=1
				ub=-ri
			else:
				ua=-ri
				ub=1
			for i in range(N):
				a,b=ab_list[i][0],ab_list[i][1]
				if (a*ua+b*ub) % pk==0:
					#print pk,n,ZZ(n[0]-n[1])
					val[i]+=ZZ(n[0]-n[1])
					#print val
	return unhit,total	
				 
def huge_test(sample_size,deg,card):
	gd=open("unhit.txt","a")
	hd=open("tested_poly.txt","a")
	for i in range(sample_size):
		q=card
		A.<t>=GF(q)['t']
		S.<x>=A['x']
		f=S.random_element(deg)
		while gcd(f.coeffs())!=1:
			f=S.random_element(deg)
		R.<X,Y>=A['X,Y']
		cl=f.coeffs()
		F=R(0)
		for j in range(f.degree()+1):
			F+=cl[j]*X^j*Y^(f.degree()-j)
		alim=10*deg		
		pairs=[[A(1),A(1)]]
		for j in range(10):
			a=A.random_element(2*deg)
			b=A.random_element(2*deg)
			if position(a)<alim and position(b)<alim:
				pairs.append([A.random_element(2*deg),A.random_element(2*deg)])
		filename="huge"+format(sample_size)+"-"+format(deg)
		makefb(f,10*deg,2*deg,filename,"human")
		u,tot=test_pair(pairs,2*deg,F,filename)
		if u!=0:
			gd.write("unhit  error "+format(f)+"\n")
		else:
			hd.write(format(f)+"\n")

	
