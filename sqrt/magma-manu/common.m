PR:=PolynomialRing;
Z:=Integers();
Q:=Rationals();
ZP<X>:=PR(Z);
x:=X;
Attach("factq.m");

// this one has wild ramification at 2.
// pol:=77*X^3-X^2+X+3;

// this one has a ramified prime of multiplicity 2 at 5, together with an
// unramified one. The ramified prime does not match a root, so can't be hit
// by a-x
// pol:=40*X^3-5*X^2+X+3;

/*
pol:=1568*x^4 - 404032074257153235*x^3 - \
10239185373556241734011564010826*x^2 +  \
97931343123749938839641205641922624916324925*x + \
32106948664745015206264092076942480719794848549316495861788;
fbsize:=1900;
idl_file := "manu/manu.ur_idl";
*/

function transform_to_ideal_valuation_vec(g,n)
	assert g[#g] eq 0;
	if n lt #g then
		assert IsZero(Vector(g[n+1..#g]));
		return Vector(g[1..n]);
	end if;
	rg:=g[1..#g-1];
	if #rg ne n then
		 rg cat:= [0 : i in [1..n-#rg]];
	end if;
	return Vector(rg);
end function;

