
// clear;

x:=Split(subdir,"/");
subdir_base:=x[#x];

idl_file := subdir cat "/" cat subdir_base cat ".ur_idl";



maxprime:=NthPrime(fbsize);

d:=Degree(pol);

// This one is the monic polynomial defining the same number field.
lcP:=LeadingCoefficient(pol);
qpol:=ZP!(Evaluate(pol, PR(Q).1/lcP)*lcP^(d-1));

disc:=Discriminant(qpol);
// the discriminant of Q is lcP^xxx * Discriminant(P)

// NF<a>:=NumberField(pol);
NF<ah>:=NumberField(qpol);
a:=ah/LeadingCoefficient(pol);

EK:=EquationOrder(NF);
// OK:=MaximalOrder(NF);

// Approximate OK the easy way.
// In any case, trial divide up to 1e6 because it's so easy.
OK:=EK;
tt:=Cputime();
printf "Finding small factors";
dfactors:=Factorization(disc :
	TrialDivisionLimit:=Maximum(10^6, maxprime),
	ECMLimit:=1000, MPQSLimit:=0);
printf "\t%o s\n", Cputime()-tt;

tt:=Cputime();
printf "Approximating MaximalOrder";
for pr in dfactors do
	p:=pr[1];
	F1:=pMaximalOrder(OK,p);
	printf " %o", Index(F1, OK);
	OK:=F1;
end for;
f:=Index(OK,EK);
printf "\nIndex: %o=%o ; %o s\n", f,Factorization(f), Cputime()-tt;

// ok, so now OK is an approximation of the maximal order, which is good enough
// to work with.
// CHEAT !
OK`Maximal:=true;

OKmat:=Matrix([Eltseq(NF!x):x in Basis(OK)]);

// We are interested in ideals of the form <p,a-q>, or more specifically in
// their cosets of index lcP which are the ideals <p, ah - lc*q>.
i:=func<p,q|ideal<OK|p,ah-LeadingCoefficient(pol)*q> >;

// Note in particular that inert primes, or more generally ideals whose
// residue class degree is >1 are discarded.

// To be on the safe side, however, we include *all* the ramified ideals, and
// all ideals whose norm divides the leading coefficient of P.
// decide whether some given ideal is of the requested form.
function is_interesting(Ip)
	if not IsPrime(Ip) then return false; end if;
	p:=Norm(Ip);
	p:=PrimeDivisors(p)[1];
	if p gt maxprime then return false; end if;
	// we start with the coarsest think we can think of.
	if IsRamified(p,OK) then return true; end if;
	// at this point, p might be a prime power.
	if lcP mod p eq 0 then return true; end if;
	if Degree(Ip) gt 1 then return false; end if;
	g:=Generators(Ip);
	assert IsPrime(p);
	assert g[1] eq OK!p;
	assert &and [g[2][i] eq 0 : i in [3..d]];
	return true;
end function;

function is_simple_ideal(Ip)
	p:=Norm(Ip);
	if not IsPrime(p) then return false; end if;
	g:=Generators(Ip);
	assert(g[1] eq OK!p);
	if #g ne 2 then return false; end if;
	return &and [g[2][i] eq 0 : i in [3..d]];
end function;

// Store all degree 1 prime ideals, plus all the ramified ones.
// In fact, it's terribly inefficient to do so (hardly 1k ideals/sec)
/*
for p in PrimesInInterval(1,maxprime) do
	try
		Fp:=Factorization(ideal<OK|p>);
		for Ip in Fp do
			if is_interesting(Ip[1]) then
				Include(~ideals, Ip[1]);
				if #ideals mod 1000 eq 0 then
					printf "%o ideals\n", #ideals;
				end if;
			end if;
		end for;
	catch e
		print "Error while handling p=", p;
		error e;	// raise error again
	end try;
end for;
*/

ideals:={@ @};
tt:=Cputime();
/* start with the special ideals */
for px in dfactors do
	p:=px[1];
	Fp:=Factorization(ideal<OK|p>);
	for Ip in Fp do
		if is_interesting(Ip[1]) then
			Include(~ideals, Ip[1]);
		end if;
	end for;
end for;

/* If we don't do so, then ordering will be screwed */
Sort(~ideals, func<a,b|Hash(a)-Hash(b)>);

special_ideals:=ideals;


tt:=Cputime()-tt;
printf "%o special ideals in %o s, %o ideals/sec\n",
	#ideals, tt, (#ideals/(tt+0.00001));

tt:=Cputime();
/* read the easy ones from precomputed C data (much much faster) */
FP:=Open(idl_file,"r");
while true do
	s:=Gets(FP);
	if IsEof(s) then
		break;
	end if;
	u:=StringToIntegerSequence(s);
	// printf "read %o -> %o \n", s, u;
	if u[1] gt maxprime then
		break;
	end if;
	I:=i(u[1], u[2]);
	if I in special_ideals or lcP mod u[1] eq 0 then
		print "Already handled special ideal: ", I;
	else
		Include(~ideals, i(u[1], u[2]));
		if #ideals mod 10000 eq 0 then
			printf "%o ideals\n", #ideals;
		end if;
	end if;
end while;
delete FP;

tt:=Cputime()-tt;
printf "%o ideals in %o s, %o ideals/sec\n",
	#ideals, tt, (#ideals/tt);


// Beware that J as an ideal in OK or J as an ideal in E are not the same
// thing.
J:=ideal<OK| [ Evaluate(pol div X^i, a) : i in [1..Degree(pol)]]>;



