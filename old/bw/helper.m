
function cf(S, i)
	return Parent(S)![Coefficient(x, i) : x in Eltseq(S)];
end function;

// submul : subtract m * [i0] to [i1]
// indices  in [0..d-1]
function smmat(d, i0, i1, m)
	r:=IdentityMatrix(Parent(M), d);
	r[1+i1][1+i0]:=-m;
	return r;
end function;

// mul [i] by m
// indices  in [0..d-1]
function mmat(d, i, m)      
	r:=IdentityMatrix(Parent(M), d);
	r[1+i][1+i]:=m;
	return r;
end function;

procedure pattern(P)
	d := Maximum([Degree(x) : x in Eltseq(P)]);
	for i in [1..Nrows(P)] do
		print "[" cat
		&cat [cf(P[i], v) eq 0 select " " else "*" : v in [0..d]]
		cat "]";
	end for;
end procedure;

