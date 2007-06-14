/**************************************************************************/
/* building relations */

if not assigned(relbatch) then
	relbatch:="";
else
	relbatch cat:= ".";
end if;

outfn:=stem cat "/" cat relbatch cat "complete";
relfn:=stem cat "/" cat relbatch cat "rel_idx";

/*
procedure handle_one(q,r,x,~rr)
	y:=r+q*x;
	fac:=Factorization(ideal<OK|a-y>);
	rel:=[];
	for I in fac do
		ix:=Index(ideals, I[1]);
		if ix eq 0 then return; end if;
		Append(~rel, <ix, I[2]>);
	end for;
	Include(~rr, rel);
	if #rr mod 100 eq 0 then
		printf "%o relations\n", #rr;
	end if;
end procedure;

FP:=Open("smallish/rel_idx","r");
// smallish/cand/achille.278125-281249.cand.rels","r");
relations:={@@};
nscanned:=0;
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
	handle_one(q,r,x,~relations);
	nscanned+:=1;
	if nscanned mod 100 eq 0 then
		printf "(%o candidates scanned)\n", nscanned;
	end if;
end while;
delete FP;
*/

/*
FR:=Open("smallish/smooth.rels","w");
for r in relations do
	fprintf FR, "%o", #r;
	for x in r do
		fprintf FR, " %o:%o", x[1]-1,x[2];
	end for;
	fprintf FR, "\n";
end for;
delete FR;
*/

	/*
FP:=Open("smallish/rel_idx","r");
IR:=Open("smallish/ab.vals","w");
rel:=1;
cand:=1;
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
	y:=r+q*x;
	idl:=ideal<OK|a-y>;
	if idl eq &* [ ideals[x[1]]^x[2] : x in relations[rel]] then
		rel +:= 1;
		fprintf IR, "%o\n", s;
	end if;
	cand +:= 1;
	if cand mod 500 eq 0 then
		printf "(%o candidates scanned, %o relations)\n", cand,rel;
	end if;
end while;
printf "(%o candidates scanned, %o relations)\n", cand,rel;
delete FP;
delete IR;
*/


rr:={};
XR:=Open(outfn,"w");
procedure handle_one(q,r,x,~ys,~nr)
	y:=r+q*x;
	if y in ys then
		return;
	end if;
	Include(~ys,y);
	fac:=Factorization(ideal<OK|a-y>);
	rel:=[];
	for I in fac do
		ix:=Index(ideals, I[1]);
		if ix eq 0 then return; end if;
		Append(~rel, <ix, I[2]>);
	end for;
	nr+:=1;
	Sort(~rel, func<a,b|a[1]-b[1]>);
	fprintf XR, "%o %o %o %o", q,r,x,#rel;
	for x in rel do
		fprintf XR, " %o:%o", x[1]-1,x[2];
	end for;
	fprintf XR, "\n";
end procedure;

tt:=Cputime();
FP:=Open(relfn,"r");
// smallish/cand/achille.278125-281249.cand.rels","r");
rel:=0;
cand:=0;
ys:={};
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
	cand+:=1;
	handle_one(q,r,x,~ys,~rel);
	if cand mod 500 eq 0 then
		printf "(%o candidates scanned, %o unique, %o relations)\n",
			cand,#ys, rel;
	end if;
end while;
Flush(XR);
delete FP;
delete XR;
tt:=Cputime()-tt;
printf "%o rels in %o s, %o rels/sec ; %o%% pass, %o cand/sec\n",
	rel, tt, (rel/(tt+0.00001)),
	Floor(100.0*rel/cand), (cand/(tt+0.00001));

// right:=i(3,0)^3*i(5,1)*i(11,6);
// left:=J*ideal<E|a-6>;
// 
// mynp:=func<p,r|[Valuation((Evaluate(Derivative(pol,k),r)/Factorial(k)),p) : k in [0..1+Degree(pol)]]>;
// 
// procedure printnps(p)
// 	for x in Roots(PR(GF(p))!pol) do
// 		printf "<%o,%o>^%o %o\n", p, x[1], x[2], mynp(p,Z!(x[1]));
// 	end for;
// end procedure;
// 
// procedure allnps()
// 	for x in TrialDivision(Discriminant(pol),10000) do
// 		printnps(x[1]);
// 	end for;
// end procedure;
// 


