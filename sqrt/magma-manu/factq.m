
intrinsic Factorization(x :: FldRatElt) -> SeqEnum, RngElt
{ the factorization of the rational x, with possibly negative multiplicities }
	fn,s:=Factorization(Numerator(x));
	fd:=Factorization(Denominator(x));
	r:=fn cat [<u[1],-u[2]>:u in fd];
	cmp:=func<a,b|Sign(a[1] - b[1])>;
	Sort(~r,cmp);
	return r, s;
end intrinsic;
intrinsic TrialDivision(x :: FldRatElt, B :: RngIntElt) -> RngIntEltFact,
	SeqEnum
{ the TrialDivision of the rational x, with possibly negative multiplicities }
	fn,cn:=TrialDivision(Numerator(x),B);
	fd,cd:=TrialDivision(Denominator(x),B);
	r:=fn cat [<u[1],-u[2]>:u in fd];
	cmp:=func<a,b|Sign(a[1] - b[1])>;
	Sort(~r,cmp);
	return r, [cn,cd];
end intrinsic;
intrinsic TrialDivision(x :: FldRatElt) -> RngIntEltFact, SeqEnum
{ the TrialDivision of the rational x, with possibly negative multiplicities }
	return TrialDivision(x,10000);
end intrinsic;
