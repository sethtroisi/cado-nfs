
// clear;

x:=Split(subdir,"/");
subdir_base:=x[#x];

// This is used in sol0.m to fetch character values ; it must end with a /
if System("[ `hostname` = 'profiterolle' ]") eq 0 then
	binpath:="main/pentium3/";
else
	binpath:="main/core2/";
end if;



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

