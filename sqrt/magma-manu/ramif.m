
_<x>:=PR(Integers());

pol:=4*x^4 - 197217352163441*x^3 - 13128578274420495258494527615369*x^2 + \
	    8689192173815010663310280676451427249558832468*x - \
		        171259963727999241986660220409340537742518155651960441457372;

function tail(P,p,r)
	/* returns P(r+p*x)=p^v*t(x) */
	coeffs:=[Integers()|];
	m:=1;
	Q:=P;
	n:=1;
	while Q ne 0 do
		Append(~coeffs,Integers()!(m*Evaluate(Q,r)));
		m*:=1/n;
		n+:=1;
		Q:=Derivative(Q);
	end while;
	// printf "%o\n", [Valuation(x,p):x in coeffs];
	v:=Minimum([Valuation(coeffs[i+1],p) + i :i in [0..#coeffs-1]]);
	// printf "coeffs: %o ; v: %o\n", coeffs, v;
	coeffs:=[Integers()| coeffs[i+1] * p^(i-v):i in [0..#coeffs-1]];
	return PR(Integers())!coeffs, v;
end function;

procedure ssets2(~l,P,p,racine,u,v,v0,kmax)
	/* We know that pol(racine+p^u * x) = p^v * P(x) */
	/* And we know that p does not divide P */
	if v ne 0 then Append(~l,<racine,p^u,v-v0>); end if;
	if v gt kmax then return; end if;
	m:=Roots(P, GF(p));
	for ra in m do
		r:=Integers()!(ra[1]);
		t,w:=tail(P,p,r);
		ssets2(~l,t,p,racine+p^u*r,u+1,v+w,v,kmax);
	end for;
end procedure;

function sievesets(P, p, kmax)
	l:=[];
	ssets2(~l,P,p,0,0,0,0,kmax);
	return l;
end function;


function lifted_root(P, Q, p, r, v)
	if v eq 1 then return r, p; end if;
	r,pv:=lifted_root(P,Q,p,r,v - (v div 2));
	pv:=pv^2 div p^(v - 2*(v div 2));
	_,m,_:=ExtendedGreatestCommonDivisor(Evaluate(Q,r),pv);
	r:=r-m*Evaluate(P,r);
	r mod:=pv;
	assert Valuation(Evaluate(P,r),p) ge v;
	return r, pv;
end function;

/* gives a description of unramified ideals above p, lifted by quadratic
 * hensel. */
function ideal_ramified(P, p, r, bound)
	kmax:=Ceiling(Log(bound)/Log(p));
	r:=lifted_root(P,Derivative(P),p,r,kmax);
	Append(~L,<p,r>);
	return L;
end function;

function ideals(P, p, bound)
	m:=Roots(P,GF(p));
	L:=[];
	kmax:=Ceiling(Log(bound)/Log(p));
	for ra in m do
		r:=Integers()!(ra[1]);
		if ra[2] gt 1 then
			i:=ideal_ramified(P,p,r,bound);
			print i;
			continue;
		end if;
		r:=lifted_root(P,Derivative(P),p,r,kmax);
		Append(~L,<p,r>);
	end for;
	return L;
end function;
