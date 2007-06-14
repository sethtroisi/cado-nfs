ye:=func<P,a,b|&+[Log(p)/p*#Roots(P,GF(p)) : p in PrimesInInterval(a,b)]>; 
sz:=func<P,a,b|[Log(AbsoluteValue(Evaluate(P,2^x)))/Log(2) : x in [a..b]]>;

function sz2(P,t)
	m:=0;
	for i in [1..100] do
		r:=Log(AbsoluteValue(Evaluate(P,Random(2^t))));
		if r gt m then m:=r; end if;
	end for;
	return m/Log(2);
end function;

nr_per_sq:=function(P,q_bits,zone_bits,i_bits,fb_bits)
	s:=sz2(P,q_bits + zone_bits + i_bits);
	y:=(ye(P,100,1000)-Log(10))/Log(2);
	s-:=y;
	u:=1.0*(s-q_bits)/fb_bits;
	prob:=Exp(u*Log(u));
	return 2*2^(zone_bits+i_bits)/prob, s;
end function;

per_sq:=function(nf,ib)
	M:=lattice_of_polys_with_root(N,root,4,2^(nf*2+ib));
	me:=pick(M);
	return nr_per_sq(me,nf,ib);
end function;

function good_poly(d,nf,ib,nt)
	// Factor base bound is expected to be approximately 2^nf
	// fbsize : #PrimesInInterval(1,2^nf)
	// ZONE macro in sieve.c : approximately 2^nf
	// -i argument to siever : approximately 2^ib
	M:=lattice_of_polys_with_root(N,root,d,2^(nf*2+ib));
	best:=0;
	for i in [1..nt] do
		me:=pick(M);
		v:=nr_per_sq(me,nf,nf,ib,nf);
		if v gt best then
			best:=v;
			best_poly:=me;
			print v;
		end if;
	end for;
	print "Size of norms (in bits)", sz2(best_poly,nf + nf + ib);
	if (LeadingCoefficient(best_poly) lt 0) then
		best_poly:=-best_poly;
	end if;
	return best_poly;
end function;


