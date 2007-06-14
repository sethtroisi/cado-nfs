
/* Lift the roots of qpol mod q^precision */
function liftroots(qpol, q, precision)
	rs1:=[Integers()|r[1] : r in Roots(qpol,GF(q))];
	rs:=[rs1];
	dqpol:=Derivative(qpol);
	i:=0;
	prc:=1;
	qq:=q;
	d:=Degree(qpol);
	while i le Log(precision)/Log(2) do
		qq:=qq^2;
		prc:=2*prc;
		i:=i+1;
		printf "Lifting roots of P mod %o^%o\n", q, prc;
		lp:=PolynomialRing(Integers(qq))!qpol;
		ldp:=PolynomialRing(Integers(qq))!dqpol;
		nrs:=rs[i];
		for j in [1..d] do
			nrs[j]:=nrs[j]-Integers()!(Evaluate(lp,nrs[j])/Evaluate(ldp, nrs[j]));
		end for;
		Append(~rs,nrs);
	end while;
	return rs;
end function;


function lift_bigprod(rs, q)
	/* Compute the big product on each residue field ; start with the
	 * highest precision, descend later on */
	i:=#rs;
	qq:=q^(2^(i-1));
	R:=Integers(qq);
	v:=[];
	for r in ChangeUniverse(rs[i],R) do
		elt_numer:=&*[(Evaluate(foo[1],r/lcP))^foo[2]:foo in numer];
		elt_denom:=&*[(Evaluate(foo[1],r/lcP))^-foo[2]:foo in denom];
		Append(~v,elt_numer/elt_denom);
	end for;
	v:=[*v*];
	while i gt 1 do
		i-:=1;
		qq:=q^(2^(i-1));
		R:=Integers(qq);
		v:=[*ChangeUniverse(v[1],R)*]cat v;
	end while;
	return v;
end function;

function roots_newton(v, q, e)
	/* Compute e-th roots by lifting -- we do it the stupid way. */
	i:=1;
	qq:=q;
	R:=Integers(qq);
	rv:=[*[R|]*];
	inv_e:=Integers()!(1/Integers(q-1)!e);
	for j in [1..d] do
		rv[1][j]:=v[1][j]^inv_e;
	end for;
	for i in [2..#v] do
		qq:=qq^2;
		printf "Lifting e-th roots mod q^%o\n",2^(i-1);
		R:=Integers(qq);
		nrv:=[R|Integers()!x: x in rv[i-1]];
		for j in [1..d] do
			r:=nrv[j];
			r +:= r/e * (1-r^e/v[i][j]);
			nrv[j]:=r;
		end for;
		Append(~rv, nrv);
	end for;
	return rv;
end function;
