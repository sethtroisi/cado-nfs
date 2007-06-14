
/**************************************************************************/
/* building solution */


gens:=[];
genfile:=stem cat "/clink.txt";
FP:=Open(genfile,"r");
while true do
	s:=Gets(FP);
	if IsEof(s) then
		break;
	end if;
	u:=StringToIntegerSequence(s);
	// printf "read %o -> %o \n", s, u;
	q:=u[1];
	r:=u[2];
	x:=u[3];
	Append(~gens, r+q*x-X);
end while;
delete FP;
printf "Read %o generating elements from %o\n", #gens, genfile;
nrel:=#gens;


W:=[];
for i in [0..block-1] do
	wfile:=Sprintf("%o/W0%o", stem, i);
	printf "Reading %o\n", wfile;
	SP:=Open(wfile, "r");
	g:=[];
	while true do
		s:=Gets(SP);
		if IsEof(s) then break; end if;
		u:=StringToIntegerSequence(s);
		Append(~g, u[1]);
	end while;
	Append(~W, g);
	delete SP;
end for;
/* also recover the result of M*W's */
MW:=[];
for i in [0..block-1] do
	wfile:=Sprintf("%o/W0%o.mul", stem, i);
	printf "Reading %o\n", wfile;
	SP:=Open(wfile, "r");
	g:=[];
	while true do
		s:=Gets(SP);
		if IsEof(s) then break; end if;
		u:=StringToIntegerSequence(s);
		Append(~g, u[1]);
	end while;
	g:=transform_to_ideal_valuation_vec(g,#ideals);
	Append(~MW, Vector(g));
	delete SP;
end for;

/*
function blob(x)
	if x gt (e div 2) then return x-e; end if;
	return x;
end function;
*/
/*
allzero:=[&and[w[i] eq 0: w in W]:i in [1..nrel]];
nz:=[[<gens[i],blob(w[i])> : i in [1..nrel] | not allzero[i]]:w in W];
*/

Attach("forge.m");

/* We need to compute elements in each of the prime ideals above e*OK ; these
 * are necessary in order to cancel the disturbing contribution of these
 * ideals when computing characters.
 */
	 /*
killers:=[];
for i in Factorization(ideal<OK|e>) do
	l:=LLL(Matrix([Eltseq(v):v in Basis(i[1])]));
	el:=Sort([OK!Eltseq(l[j]):j in [1..Nrows(l)]],func<a,b|Norm(a)-Norm(b)>);
	Append(~killers,<Index(ideals,i[1]),el[1]>);
end for;
*/

print "Computing characters";

cfilename := stem cat "/" cat stem cat ".characters";

/* First decide which characters are usable */

/* schirokauer-type characters */


pAdicPrecision:=10;

charstring:="";
if not IsRamified(e,OK) and #Roots(PR(GF(e))!pol) eq 0 then
	flm:=Lcm([Degree(i[1]):i in Factorization(ideal<OK|e>)]);
	charstring cat:= Sprintf("sf%o",flm);
else for fac in Factorization(PR(pAdicRing(e,pAdicPrecision))!pol) do
	/* Both in the e>1 and f>1 cases, we are lazy here ; clearly the C
	 * program should be augmented for that. */
	if Degree(fac[1]) gt 1 then continue; end if;
	r:=Roots(fac[1]);
	assert #r eq 1 and r[1][2] eq 1;
	r:=r[1][1];
	charstring cat:= Sprintf(" sr%o",Z!r);
end for;
end if;


Qs:=[];
Qc:=[];
q:=1;
while #Qs lt 20 do
	repeat q+:=e; until IsPrime(q);
	for u in Roots(PR(pAdicRing(q,pAdicPrecision))!pol) do
		Append(~Qs, Sprintf("%or%o", q, Z!u[1]));
		Append(~Qc, <q,Z!u[1]>);
	end for;
end while;
// for v in Qs do charstring cat:= " " cat v; end for;



/* This does not work as exepected. Weird */
	/*
function charfile()
	try
		return Open(cfilename,"r");
	catch err
		printf "Character file %o does not exist\n", cfilename;
		print "Generating";
		polyfilename:=stem cat "/" cat stem cat ".poly";
		cmd:="make " cat polyfilename;
		print cmd; System(cmd);
		flm:=Lcm([Degree(i[1]):i in Factorization(ideal<OK|e>)]);
		cmd:=Sprintf("core2/characters %o/clink.txt %o %o < %o > %o",
			stem, e, flm, polyfilename, cfilename);
		print cmd; System(cmd);
		System("sleep 1");
	end try;
	print "Hello, world";
	CF:=Open(cfilename,"r");
	CF;
	return CF;
end function;
CF:=charfile();
*/

if System("[ `hostname` = 'profiterolle' ]") eq 0 then
	machinename:="pentium3";
else
	machinename:="core2";
end if;

polyfilename:=stem cat "/" cat stem cat ".poly";
cmd:="make " cat polyfilename;
print cmd; System(cmd);
cmd:=Sprintf("if [ ! -f %o ] ; then %o/characters %o/clink.txt %o %o < %o > %o ; fi",
	cfilename, machinename, stem, e, charstring, polyfilename, cfilename);
print cmd; System(cmd);

CF:=Open(cfilename,"r");
print "Reading characters file";
genchars:=[];
while true do
	s:=Gets(CF);
	if IsEof(s) then break; end if;
	u:=StringToIntegerSequence(s);
	Append(~genchars, Vector(GF(e),u));
end while;
delete CF;


/*
genchars:=[additive_characters(Evaluate(g,a),e):g in gens];
*/
 
// prods:=[&*[Evaluate(gens[i],a)^(Z!w[i]): i in [1..nrel]]  : w in W];


charblock:=[ &+ [ w[i]*genchars[i] : i in [1..nrel] ] : w in W];

ns:=Matrix(Basis(Nullspace(Matrix(charblock))));
nsz:=Matrix(Nrows(ns),Ncols(ns),ChangeUniverse(Eltseq(ns),Integers()));
nsLLL:=LLL(VerticalJoin(nsz, e*IdentityMatrix(IntegerRing(),Ncols(ns))));

xdeps:=[];
for i in [1..Nrows(nsLLL)] do
	if nsLLL[i] ne 0 then
		Append(~xdeps, nsLLL[i]);
	end if;
end for;
xdeps:=Sort(xdeps, func<a,b|Norm(a)-Norm(b)>);

print "Dependencies", xdeps;

/* Now try to obtain the decomposition of our p-th power */

depnum:=0;
decomp:=[];
while #[u:u in decomp|u[2] ne 0] eq 0 do
	depnum +:=1;
	if depnum gt #xdeps then break; end if;
	chosen_chardep:=xdeps[depnum];
	decomp:=[<gens[i],&+[m[j]*W[j][i] : j in [1..#W]]> : i in [1..nrel]]
		where m is chosen_chardep;
end while;
if depnum gt #xdeps then
	print "Argh, only trivial dependencies found\n";
end if;
printf "Chosen dependency #%o : %o\n", depnum, chosen_chardep;



function ispower_quick_product(pform,e)
	assert not IsRamified(e,OK);
	qs:=[];
	q:=1+e*10^6;
	while #qs lt 1 do
		if IsPrime(q) and IsInert(q,OK) /* not IsRamified(q,OK) */ then
			Append(~qs,q);
		end if;
		q+:=e;
	end while;
	vals:=[* GF(q)!1 : q in qs *];
	for i in [1..#qs] do
		q:=qs[i];
		// fm:=Lcm([Degree(f[1]) : f in
		// Factorization(PolynomialRing(GF(q))!pol)]);
		v:=GF(q^d)!1;
		ex:=(q^d-1)div e;
		for j in [1..#pform] do
			x:=pform[j];
			y:=NF!x[1];
			rn:=power_element_mod(mod_out(y, q),ex,q);
			z:=GF(q^d)!ChangeUniverse(Eltseq(rn),GF(q));
			v *:= z^x[2];
			if j mod 100 eq 0 then print j; end if;
		end for;
		vals[i]:=v;
	end for;
	return &and[IsOne(x):x in vals], vals;
end function;

function ispower_quick_quotient(num,den,e)
	assert not IsRamified(e,OK);
	qs:=[];
	q:=1+e*10^6;
	while #qs lt 10 do
		if IsPrime(q) and not IsRamified(q,OK) then
			Append(~qs,q);
		end if;
		q+:=e;
	end while;
	for q in qs do
		fm:=Lcm([Degree(f[1]) : f in
		Factorization(PolynomialRing(GF(q))!pol)]);
		rn:=power_element_mod(mod_out(NF!num, q),(q^fm-1)div e,q);
		rd:=power_element_mod(mod_out(NF!den, q),(q^fm-1)div e,q);
		if rn ne rd then return false,q,mod_out(rn/rd,q); end if;
	end for;
	return true;
end function;
function ispower_quick(elt,e)
	return ispower_quick_quotient(elt,Parent(elt)!1,e);
end function;

// Q:=[]; q:=10^9; while #Q lt 10 do q:=NextPrime(q); Append(~Q,q); end while;


print "Computing logarithmic embeddings of the root";
/* This could be made with mpfr instead */

RR:=RealField(real_precision);
CC:=ComplexField(real_precision);

sqrt_2_RR:=Sqrt(RR!2);

embspace:=VectorSpace(RR,d);

rootsCC:=[x[1]:x in Roots(qpol, CC)];

conjsign:= [ IsReal(x) select 0 else Imaginary(x) gt 0 select 1 else -1
	: x in rootsCC ];

function MyConjugates(b)
	return [CC|Evaluate(PR(Q)!Eltseq(NF!b),u):u in rootsCC];
end function;

/* Prepare the matrix for bounding coefficients in terms of logarithmic
 * embeddings ; the coefficients we look for are those in the ah^i polynomial
 * basis.
 */

cmat:=Matrix([Eltseq(Evaluate(pol,T) div ((T-aj)*Evaluate(Derivative(qpol),aj))):aj in rootsCC]) where T is PolynomialRing(CC).1;

function bound_from_logembs(v)
	return [&+[AbsoluteValue(cmat[i][j])*Exp(v[j]):j in [1..d]]:i in [1..d]];
end function;


function logembeddings(v)
	return embspace![Log(AbsoluteValue(w)):w in MyConjugates(v)];
end function;

embs0:=1/e * &+[embspace| x[2] * logembeddings(Evaluate(x[1],a)):x in decomp];

/* These might be used later */
numer:=[x:x in decomp|x[2] gt 0];
denom:=[x:x in decomp|x[2] lt 0];


	/* computes the ideal valuation vector */
function ivec(elt)
	VS:=RSpace(Integers(),#ideals);
	v:=VS!0;
	I:=ideal<OK|elt>;
	f:=Factorization(I);
	for i in f do
		idx:=Index(ideals,i[1]);
		assert idx gt 0;
		v[idx]:=i[2];
	end for;
	return v;
end function;

function ivec0(elt)
	I:=ideal<OK|elt>;
	return Vector([Valuation(I,ideals[i]):i in [1..10]]);
end function;

print "Computing ideal valuations";

vals:=&+[m[j]*MW[j] : j in [1..#W]] where m is chosen_chardep;

// this is the version that does not rely on the C files.
// vals:=&+[ d[2]*ivec(Evaluate(d[1],a)) : d in decomp];
assert IsZero(VectorSpace(GF(e),#ideals)!vals);
vals:=[v div e: v in Eltseq(vals)];

/* vals_plus is vals[1], vals_minus is vals[2] */
vals0:=[[v gt 0 select v else 0:v in vals],[v lt 0 select -v else 0:v in vals]];


/* Summarize the norm contribution from ideals, separately in the numerator
 * and denominator */

for i in ideals do assert(Norm(i) gt 1); end for;
lognorms0:=[&+[w[i]*Log(RR!Norm(ideals[i])): i in [1..#ideals]]:w in vals0];

/*
lognorm_numer:=&+[RR| x[2]*Log(AbsoluteValue(Norm(Evaluate(x[1],a)))):x in decomp|x[2] gt 0];
lognorm_denom:=&+[RR|-x[2]*Log(AbsoluteValue(Norm(Evaluate(x[1],a)))):x in decomp|x[2] lt 0];
emb_numer:=&+[VectorSpace(RR,d)|x[2]*VectorSpace(RR,d)![Log(AbsoluteValue(v)) : v in MyConjugates(Evaluate(x[1],a))]:x in numer];
emb_denom:=-&+[VectorSpace(RR,d)|x[2]*VectorSpace(RR,d)![Log(AbsoluteValue(v)) : v in MyConjugates(Evaluate(x[1],a))]:x in denom];
emb_quotient:=emb_numer - emb_denom;
*/


