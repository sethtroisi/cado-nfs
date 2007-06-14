
intrinsic Lfunc(N :: RngIntElt,a :: FldRatElt,c :: FldReElt) -> FldReElt
{ Pomerance's L function }
	return Exp(c*Log(N)^a*Log(Log(N))^(1-a));
end intrinsic;

intrinsic give_data(N :: RngIntElt,d :: RngIntElt,B :: RngIntElt,X :: RngIntElt,skew :: RngIntElt)
{ prints some data }
	if skew gt X then
		print "[ this skew value is too large ]";
	end if;
	norms1:=N^(1/(d+1)) * Maximum(skew,X)^(d/2);
	prob1:=u^u where u is (Log(norms1)/Log(B));
	printf "FB bound: %o [%o bits]\n", Floor(B), Ilog2(B);
	printf "sieve bound: %o [%o bits]\n", Floor(X), Ilog2(X);
	printf "matrix size: %o\n", Floor(B/Log(B));
	printf "max skew: %o bits\n", Ilog2(N^(1/(d+1)));
	printf "(here) size of norms: %o bits\n", Ilog2(norms1);
	printf "(here) smoothness prob: 1/%o\n", Floor(prob1);
	printf "estimated #rels obtained: %o (%o%%)\n",
		Floor(2*X/prob1), Floor(2*X/prob1/(B/Log(B))*100);
end intrinsic;

intrinsic hints(N :: RngIntElt)
{ hints }
	printf "Hints for %o [%o bits]\n", N, Ilog2(N);
	c0:=(3/2)^(1/3);
	c1:=(32/9)^(1/3);
	c2:=(4/9)^(1/3);

	d:=c0*(Log(N)/Log(Log(N)))^(1/3) - 1;
	printf "degree: %.3o --> %o\n", d, Round(d);
	d:=Round(d);

	X:=Floor(Lfunc(N, 1/3, c1));
	B:=Floor(Lfunc(N, 1/3, c2));
	give_data(N,d,B,X,X);

	norms2:=Lfunc(N,2/3,1/c0+c1*c0/2);
	prob2:=Lfunc(N, 1/3, 1/3*(1/c0+c1*c0/2)/c2);
	printf "(theory) size of norms: %o bits\n", Ilog2(norms2);
	printf "(theory) smoothness prob: 1/%o\n", Floor(prob2);
end intrinsic;

intrinsic score(P :: RngUPolElt) -> FldReElt
{}
	return &+ [ #Roots(P, GF(p)) * Log(p)/p : p in PrimesInInterval(1,100)
		];
end intrinsic;

intrinsic padvect(p :: RngUPolElt, d :: RngIntElt) -> SeqEnum
{ same as EltSeq, but with specified length }
	v:=Eltseq(p);
	if #v lt d+1 then
		v cat:= [0: i in [#v..d]];
	end if;
	return v;
end intrinsic;

intrinsic padvect(p :: RngIntElt, d :: RngIntElt) -> SeqEnum
{ same as EltSeq, but with specified length }
	return padvect(PolynomialRing(Integers())!p, d);
end intrinsic;

intrinsic scale_up(~M :: MtrxSpcElt, s :: FldRatElt)
{ }
	a:=AbsoluteValue(1/s);
	if a gt 1 then
		a:=Floor(a);
		for i in [1..Nrows(M)] do
		for j in [1..Ncols(M)] do
			M[i][j] *:= a^(Ncols(M)-j);
		end for;
		end for;
	else
		a:=Floor(1/a);
		for i in [1..Nrows(M)] do
		for j in [1..Ncols(M)] do
			M[i][j] *:= a^(j-1);
		end for;
		end for;
	end if;
end intrinsic;
intrinsic scale_down(~M :: MtrxSpcElt, s :: RngIntElt)
{ }
	scale_down(~M, s/1);
end intrinsic;

intrinsic scale_down(~M :: MtrxSpcElt, s :: FldRatElt)
{ }
	a:=AbsoluteValue(1/s);
	if a gt 1 then
		a:=Floor(a);
		for i in [1..Nrows(M)] do
		for j in [1..Ncols(M)] do
			assert M[i][j] mod a^(Ncols(M)-j) eq 0;
			M[i][j] div:= a^(Ncols(M)-j);
		end for;
		end for;
	else
		a:=Floor(1/a);
		for i in [1..Nrows(M)] do
		for j in [1..Ncols(M)] do
			assert M[i][j] mod a^(j-1) eq 0;
			M[i][j] div:= a^(j-1);
		end for;
		end for;
	end if;
end intrinsic;
intrinsic scale_up(~M :: MtrxSpcElt, s :: RngIntElt)
{ }
	scale_up(~M, s/1);
end intrinsic;

intrinsic lattice_of_polys_with_root(N :: RngIntElt, root :: RngIntElt, d :: RngIntElt, skew :: RngIntElt) -> MtrxSpcElt
{ returns the lattice of admissible polys }
	vlist:=[];
	Append(~vlist, padvect(N,d));
	X:=PolynomialRing(Integers()).1;
	for i in [1..d] do
		Append(~vlist, padvect(X^i - root^i mod N, d));
	end for;

	Mv := Matrix(vlist);
	// print Mv;
	scale_up(~Mv, skew);
	Mv2:=LLL(Mv);
	scale_down(~Mv2, skew);
	return Mv2;
end intrinsic;

intrinsic pick(M :: MtrxSpcElt) -> RngUPolElt
{}
	ZP:=PolynomialRing(Integers());
	seen:={@@};
	w:=2;
	repeat
		foo:=Vector([Random([-w..w]) : i in [1..Nrows(M)]]);
		if foo in seen then
			loop +:=1;
			if loop ge 100 then
				print "looping, w=", w, "seen=", #seen;
				w+:=1;
				loop := 0;
			end if;
			continue;
		end if;
		loop := 0;
		Include(~seen, foo);
		P:=ZP!Eltseq(foo*M);
	until Degree(P) eq Nrows(M)-1 and
		/* IsMonic(P) and */
			/*
		(IsPrime(x) or x eq 1 where x is LeadingCoefficient(P)) and
			*/
			Maximum([1] cat PrimeDivisors(LeadingCoefficient(P))) le 100
				and
			IsIrreducible(P);
	//print "real roots :", Roots(P, RealField());
	return P;
end intrinsic;

intrinsic pick_poly(N :: RngIntElt, root :: RngIntElt, d :: RngIntElt, skew :: RngIntElt) -> RngUPolElt, MtrxSpcElt
{ picks a polynomial }
	Mv2:=lattice_of_polys_with_root(N, root, d, skew);
	P:=pick(Mv2);
	assert Evaluate(P, root) mod N eq 0;
	return P, Mv2;
end intrinsic;

NFSFORGE:=recformat<
	N : RngIntElt,
	m1 : RngIntElt,
	m2 : RngIntElt,
	root : RngIntElt,
	
	d : RngIntElt,
	Xbound : RngIntElt,
	B : RngIntElt,

	P : RngUPolElt,
	I : SeqEnum,
	nothanks : RngIntElt,

	skew : RngIntElt,

	ncharacters : RngIntElt,

	lines : SetIndx
>;
	
intrinsic describe(NFS :: Rec)
{ prints sort of a job file }
	assert Names(NFS) eq Names(NFSFORGE);
	printf "n: %o\n", NFS`N;
	printf "m: \n";
	for i in [NFS`d .. 0 by -1] do
		printf "c%o: %o\n", i, Coefficient(NFS`P, i);
	end for;
	printf "Y0: %o\n", NFS`m2;
	printf "Y1: %o\n", -NFS`m1;
	printf "skew: %o\n", NFS`skew;
	printf "rlim: %o\n", NFS`B;
	printf "alim: %o\n", NFS`B;
	printf "lpbr: 0\n";
	printf "lpba: 0\n";
	printf "mpbr: 0\n";
	printf "mpba: 0\n";
	printf "rlambda: 2.5\n";
	printf "alambda: 2.5\n";
	printf "a0: %o\n", -NFS`Xbound;
	printf "a1: %o\n", NFS`Xbound;
	printf "b0: %o\n", 1;
	printf "b1: %o\n", 1;
end intrinsic;


intrinsic SmoothNorm(NFS :: Rec, q :: RngIntElt) -> BoolElt, SeqEnum
{ tests if alpha+q gives a relation }
	assert Names(NFS) eq Names(NFSFORGE);
	// principal ideal generated by alpha+q.
	n:=Evaluate(NFS`P, -q);
	f:=Factorization(n);
	if not &and [ x[1] le NFS`B : x in f ] then return false,_; end if;
	if not &and [ NFS`nothanks mod x[1] ne 0 : x in f ] then return false,_; end if;
	fl:=[];
	for i in [1..#NFS`I] do
		p:=NFS`I[i][1];
		rp:=NFS`I[i][2];
		v:=Valuation(n, p);
		if v eq 0 then continue; end if;
		if GF(p)!-q eq rp then
			Append(~fl, [ i, v]);
			assert n mod p^v eq 0;
			n div:= p^v;
		end if;
	end for;
	if AbsoluteValue(n) ne 1 then
		print "argh", q, n;
	end if;
	return true, [#fl] cat &cat fl;
end intrinsic;

intrinsic try_it(~NFS :: Rec, q :: RngIntElt)
{ one trial }
	assert Names(NFS) eq Names(NFSFORGE);
	x,xf:=SmoothNorm(NFS, q);
	if x and <q,xf> notin NFS`lines then
		Include(~NFS`lines, <q,xf>);
		print Sprintf("%o/%o", #NFS`lines, #NFS`I), q, xf;
	end if;
end intrinsic;

intrinsic try_from(~NFS :: Rec,x0 ::RngIntElt)
{ everything from x0 }
	assert Names(NFS) eq Names(NFSFORGE);
	q:=x0;
	while q le 2*NFS`Xbound do
		if (q mod 1000 eq 0) then
			printf "tested %o/%o\n", q, (2*NFS`Xbound);
		end if;
		try_it(~NFS, (q div 2) - (q mod 2) * q);
		if #NFS`lines gt #NFS`I+1+NFS`ncharacters then break; end if;
		q+:=1;
	end while;
end intrinsic;

intrinsic mod_out(x ::FldNumElt,p ::RngIntElt) -> FldNumElt
{  mods out the coeffs }
        return Parent(x)![ Integers()!Integers(p)!u : u in Eltseq(x) ] ;
end intrinsic;

intrinsic power_element_mod(x :: FldNumElt,n :: RngIntElt,p :: RngIntElt) -> FldNumElt
{ powers an element in a number field, modding out the coefficients by p }
        if n eq 0 then return Parent(x)!1; end if;
        y := mod_out(power_element_mod(x, n div 2, p),p);
        y := mod_out(y*y,p);
        if n mod 2 eq 0 then return y; else return mod_out(x*y,p); end if;
end intrinsic;
intrinsic bexpand(n :: RngIntElt,p :: RngIntElt,k :: RngIntElt) -> SeqEnum
{ expand n in base p, gives k coeffs }
        L:=[Integers()|];
        for i in [0..k-1] do
                Append(~L,n mod p);
                n div:=p;
        end for;
        assert n eq 0;
        return L;
end intrinsic;

intrinsic additive_characters(x :: FldNumElt,p :: RngIntElt) -> SeqEnum
{ schirokauer maps }
        NF:=Parent(x);
        P:=DefiningPolynomial(NF);
        assert(IsPrime(p));
        m:=Lcm([Degree(f[1]) : f in Factorization(PolynomialRing(GF(p))!P)]);
  	k:=1;	/* k>1 is wrong */
        g:=power_element_mod(x,p^m-1,p^(k+1))-1;
	VS:=VectorSpace(GF(p),Degree(NF));
	return VS![GF(p)|(Integers()!u div p): u in Eltseq(g) ];
	/*
	r:=VS!0;
        s:=[ bexpand(Integers()!u div p,p,k) : u in Eltseq(g) ];
        return s;
  */
end intrinsic;

intrinsic additive_characters(xx::SeqEnum,p :: RngIntElt) -> SeqEnum
{ schirokauer maps }
	NF:=Parent(xx[1][1]);
	P:=DefiningPolynomial(NF);
	assert(IsPrime(p));
	m:=Lcm([Degree(f[1]) : f in Factorization(PolynomialRing(GF(p))!P)]);
	VS:=VectorSpace(GF(p),Degree(NF));
	r:=VS!0;
	for x in xx do
		g:=power_element_mod(x[1],p^m-1,p^2)-1;
		r+:=x[2]*VS![GF(p)|(Integers()!u div p): u in Eltseq(g) ];
	end for;
	return Eltseq(r);
end intrinsic;

intrinsic additive_characters_1(x :: FldNumElt,p :: RngIntElt,k :: RngIntElt) -> SeqEnum
{ schirokauer maps }
	require k eq 1 : "character lifting not well defined";
	return additive_characters(x,p);
end intrinsic;

intrinsic fixup(v,e) -> SeqEnum
{}
	Z:=Integers();
	return [Z|(e-x lt x select -(e-x) else x) where x is Z!x0 : x0 in
		Eltseq(v)];
end intrinsic;
