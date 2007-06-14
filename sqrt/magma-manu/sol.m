load "common.m";

load "load_medium38.m";
e:=23;
block:=8;
real_precision:=1000;

load "setup.m";
load "sol0.m";
load "icheck.m";
load "sol_init.m";
load "onestep.m";
load "lifting.m";
procedure info(nsteps, ~c)
	delta:=(c`lognorms[1]-c`lognorms[2])/d;
	now:=[Round(x):x in c`lognorms]
		cat [Round(x-delta) : x in Eltseq(c`embs)];
	printf "Step %o : %o\n", nsteps, now;
end procedure;
tt:=Cputime();
c:=c0;
s:=2;
nsteps:=0;
info(nsteps, ~c);
while not assigned(c`done) do
	s:=3-s;onestep(s,~c);
	nsteps+:=1;
	info(nsteps, ~c);
end while;
tt:=Cputime()-tt;
printf "Done %o steps, %o seconds\n", nsteps, tt;


/* This should be the ``denominator'' of OK ; smallest d such that dOK is
 * contained  in the equation order */
OKmat:=Matrix([Eltseq(NF!x):x in Basis(OK)]);
multiplier:=Lcm([Denominator(x):x in Eltseq(OKmat)]);

cbounds:=bound_from_logembs(Eltseq(c`embs));
cbmax:=Maximum(cbounds);
print "bound on largest coeff of big/root^e : ", cbmax;
/* Take some security margin */
cbmax*:=multiplier;
cbmax*:=2000;



/* This would be the approach where we select a set of several primes to
 * compute via CRT.
Q:=[]; q:=10^9;
while #Q lt Ceiling(Log(cbmax)/9) do
	q:=NextPrime(q);
	Append(~Q,q);
end while;
 */

/* Since we prefer to work with Hensel lifting instead, we choose one single
 * prime ; to ease things a bit, we choose a completely split one */

q:=10^9; 
repeat q:=NextPrime(q);
until #Roots(pol,GF(q)) eq d and Gcd(e,q-1) eq 1;

precision:=Ceiling(Log(cbmax)/Log(q));
/* Exaggerate precision */
precision+:=4;
precision*:=2;

rs:=liftroots(qpol, q, precision);
v:=lift_bigprod(rs, q);
rv:=roots_newton(v, q, e);

/*
maps:=[[hom<NF->Integers(q^(2^(i-1)))|r>:r in rs[i]]:i in [1..#rs]];
[v[i] eq [m(big_product):m in maps[i]] : i in [1..#rs]];
[rv[i] eq [m(totalroot):m in maps[i]] : i in [1..#rs]];
*/

/* From now on, we focus only on the highest precision approximations */

qq:=q^precision;
R:=Integers(qq);
/* maybe precision is not a power of 2 */
high_roots:=[R|Integers()!x : x in rv[#rv]];
high_ahs:=[R|Integers()!x : x in rs[#rs]];

proddeltas_modqq:=[
	&*[Parent(r)|Evaluate(PolynomialRing(R)!Eltseq(NF!(d[1])),r)^d[2]:d in c`deltas]
		:r in high_ahs];

/*
proddeltas_modqq eq [R|&*[m(d[1])^d[2]:d in c`deltas]:m in maps[#maps]];
[high_roots[j]/proddeltas_modqq[j]:j in [1..d]] eq [R|m(rtail):m in maps[#maps]];
*/

/* Now the images of the mapping are computed by a transformation by a
 * vandermonde matrix; recover the real thing now ! */

images:=Vector([R|high_roots[j]/proddeltas_modqq[j]:j in [1..d]]);
vdm:=Transpose(Matrix([[R|r^j:j in [0..d-1]]:r in rs[#rs]]));

coeffs:=ChangeUniverse(Eltseq(multiplier*images*vdm^-1),Integers());

for i in [1..#coeffs] do
	cf:=coeffs[i];
	if cf gt qq-qq/1000 then
		cf-:=qq;
	end if;
	if cf gt qq/1000 then
		print "Argh, lifting error";
	end if;
	coeffs[i]:=cf;
end for;

root_in_NF:=NF!coeffs/multiplier;

/* Now compute the embedding mod N */

assert Evaluate(pol,P) mod N eq 0;

mapN:=hom<NF->Integers(N)|P*lcP>;

// product_modN:=&*[(Integers(N)!Evaluate(d[1],P))^d[2] : d in decomp];
product_modN:=&*[mapN(Evaluate(d[1],a))^d[2]:d in decomp];
product_modN;

rootmodN:=&*[mapN(d[1])^d[2]:d in c`deltas] * mapN(root_in_NF);

tg,mtg:=TorsionUnitGroup(NF);
for v in tg do
	u:=mtg(v);
	print mapN(u)*rootmodN^e;
end for;


procedure ucheck(c,root_in_NF)
	NF:=Parent(root_in_NF);
	assert Evaluate(pol,P) mod N eq 0;
	mapN:=hom<NF->Integers(N)|P*lcP>;

	// product_modN:=&*[(Integers(N)!Evaluate(d[1],P))^d[2] : d in decomp];
	product_modN:=&*[mapN(Evaluate(d[1],a))^d[2]:d in decomp];
	product_modN;

	rootmodN:=&*[mapN(d[1])^d[2]:d in c`deltas] * mapN(root_in_NF);

	tg,mtg:=TorsionUnitGroup(NF);
	for v in tg do
		u:=mtg(v);
		print mapN(u)*rootmodN^e;
	end for;
end procedure;


/*
elt_numer:=&*[(foo[1]-a)^foo[2]:foo in numer];
elt_denom:=&*[(foo[1]-a)^-foo[2]:foo in denom];
foo:=elt_numer/elt_denom * (&*[x[1]^x[2]:x in c`deltas])^e;

cc:=MyConjugates(foo);
[Floor(Log(AbsoluteValue(x))/e):x in cc];
v[1]-v[2] where v is Matrix(c`embs);

nf:=Norm(foo);
Log(AbsoluteValue(nf))/e;
c`lognorm[1] - c`lognorm[2];
c`lognorm;

multipliers:=[i: i in [10^6..10^6+100]|IsPrime(i*e+1)];
q:=multipliers[1]*e+1;
fm:=Lcm([Degree(f[1]) : f in Factorization(PolynomialRing(GF(q))!pol)]);
power_element_mod(mod_out(foo, q),(q^fm-1)div e,q);



*/

/*
M2:=[];
c:=(Lbound/fullnorm*Discriminant(pol/lcP))^(1/d);
for i in [1..Nrows(M)] do
	v:=Eltseq(M[i]);
	embs:=VectorSpace(RR,d)![Log(AbsoluteValue(w)) : w in MyConjugates(OK!v)];
	embs-:=emb_quotient/e;
	v cat:= [ Floor(c*x) : x in Eltseq(embs) ];
	Append(~M2, v);
end for;

v:=M[1];
lengthv:=Sqrt(Norm(Vector(v)));
Norm(OK!Eltseq(v))/fullnorm le C;
*/


/*
lengthv le 2^((d-1)/4)*fullnorm^(1/d);
Determinant(SylvesterMatrix(ZP!Eltseq(v),pol)) eq Norm(OK!Eltseq(v));
*/


/*
> foo:=Matrix(5,5,[89, 363, 476,0,0,0,89, 363, 476,0,0,0,89, 363, 476,10085875\
669933844, -12581047, 1, 1,0,0,10085875669933844, -12581047, 1, 1]);
> Determinant(foo);
10971047087983548275457460933104725151499
> $1/lengthv^d*lengthf^(d-1);
5.03477875830162828896328340352E63
> Determinant(foo);           
10971047087983548275457460933104725151499
> $1/(lengthv^d*lengthf^(d-1));
0.486548289687611068194098111924
	*/
	



/*
elt_numer:=&*[(foo[1]-a)^foo[2]:foo in numer];
elt_denom:=&*[(foo[1]-a)^-foo[2]:foo in denom];

In:=ideal<OK|elt_numer>;
Id:=ideal<OK|elt_denom>;          
[Valuation(In,i) : i in ideals[1..100]];
[Valuation(Id,i) : i in ideals[1..100]];


biggens:=[&*[(foo[1]-a)^foo[2]:foo in z0] : z0 in nz];
really_power:=&*[biggens[i]^foo[i] : i in [1..Ncols(ns)]] where foo is xdeps[1];
*/
/*
while c`lvecs[1][3] ge 30 and c`lvecs[2][3] ge 30 do
onestep(1,~c);c`lognorm;[[Floor(y):y in Eltseq(x)]:x in c`embs];
onestep(2,~c);c`lognorm;[[Floor(y):y in Eltseq(x)]:x in c`embs];
end while;
*/


/*
elt_numer:=&*[(Evaluate(foo[1],a))^foo[2]:foo in numer];
elt_denom:=&*[(Evaluate(foo[1],a))^-foo[2]:foo in denom];
eltx:=[elt_numer,elt_denom];
eltx0:=eltx;
big_product:=elt_numer/elt_denom;
*/


/*
	// expensive version.
	c:=c0;
	c`I:=[
		&*[ideals[i]^vals[i]:i in [1..#ideals]|vals[i] gt 0],
		&*[ideals[i]^-vals[i]:i in [1..#ideals]|vals[i] lt 0]];
	c1:=c;

	elt_numer:=&*[(lcP*Evaluate(foo[1],a))^foo[2]:foo in numer];
	elt_denom:=&*[(lcP*Evaluate(foo[1],a))^-foo[2]:foo in denom];
	assert elt_numer in OK;
	assert elt_denom in OK;
	xp:=&+[foo[2]:foo in numer];
	xm:=&+[-foo[2]:foo in denom];
	if xp gt xm then
		elt_denom *:= lcP^(xp-xm);
	else
		elt_numer *:= lcP^(xm-xp);
	end if;
	// These (elt_numer, elt_denom) are plain algebraic integers

	big_product:=elt_numer/elt_denom;

	s:=2;
	c:=c0;
	nsteps:=0;
	tt:=Cputime();
	tail:=big_product;

	I3s:=ideals[1..3];

	for i in [1..43] do
	s:=3-s;onestep(s,~c);
	tail/:=(s[1]^(e*s[2])) where s is c`deltas[#c`deltas];
	Vector([Valuation(d[1],u):u in I3s]) where d is c`deltas[#c`deltas];
	&+[RSpace(Z,3)|Vector([d[2]*Valuation(d[1],u):u in I3s]):d in c`deltas];
	1/e*Vector([Rationals()|Valuation(vvv,u):u in I3s]) where vvv is ideal<OK|tail>;
	end for;




	// computed root:
	croot:=NF!coeffs * &*[x[1]^x[2]:x in c`deltas];

//	NF!coeffs eq rtail;


	// Descend manually to zero ; or figure out a bound on the
	// coefficients of the tail given the logarithmic embeddings
	root_numer:=&*[x[1]:x in c`deltas | x[2] gt 0];
	root_denom:=&*[x[1]:x in c`deltas | x[2] lt 0];
	rootx:=[root_numer,root_denom];
	root:=&*[x[1]^x[2]:x in c`deltas];
	re:=root^e;
	tail:=big_product/re;
	ttt,rtail:=IsPower(tail,e);
	totalroot:=root*rtail;

*/
