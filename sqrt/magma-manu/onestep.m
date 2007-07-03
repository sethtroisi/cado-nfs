
sign:=[1,-1];
s:=1;

procedure onestep(s, ~c)
	sign:=[1,-1];

	if assigned(c`done) then
		print "finished";
		return;
	end if;

	/* These two must be close to eachother ! */
		/*
	print [RealField(10)|&+Eltseq(c`embs), c`lognorms[1]-c`lognorms[2]];
	*/

	selected, fullnorm:=pick_ideals(s,c);
	I:=c`H[s]*&*{Parent(c`H[s])|ideals[foo[1]]^foo[2]:foo in selected};
	M:=LLL(Matrix([Eltseq(x):x in Basis(I)]));

	/* Use det(M) instead of Lbound -- just for a try */
	// cst  := Log(AbsoluteValue(Determinant(M)));
	cst  := Log(RR!Lbound);
	cst -:= Log(RR!fullnorm);

	// This term is not there I think.
	// cst -:= (d-1)*Log(AbsoluteValue(lcP));

	// This one really is the discriminant of the ORDER !
	cst -:= 1/2 * Log(RR!AbsoluteValue(Discriminant(OK)));
	// cst -:= 1/2 * Log(AbsoluteValue(Discriminant(qpol)));
	
	// could also write:
	// Determinant(OKmat)*Sqrt(AbsoluteValue(Discriminant(qpol)))
	// but it's the same (the determinant being with respect to the basis
	// of the equation order.


	cst +:= (c`lognorms[1]-c`lognorms[2]) * sign[s];
	cst /:= d;

	/* equation middle of page 10. */
	assert Sqrt(Norm(Vector(M[1]))) le 2^((d-1)/4)*fullnorm^(1/d);

	
 	M2:=[];
 	scales:=[Exp(cst-sign[s]*c`embs[i]) : i in [1..d]];
 	for i in [1..Nrows(M)] do
 		v:=Eltseq(M[i]);
 		// embs:=[Log(AbsoluteValue(w)) : w in MyConjugates(OK!v)];
 		embs:=MyConjugates(OK!v);
 		for i in [1..d] do
 			if conjsign[i] eq 0 then
 				continue;
 			elif conjsign[i] eq 1 then
 				embs[i]:=sqrt_2_RR*Real(embs[i]);
 			elif conjsign[i] eq -1 then
 				embs[i]:=sqrt_2_RR*Imaginary(embs[i]);
 			end if;
 		end for;
 		ChangeUniverse(~embs, RR);
 		v cat:= Eltseq(embs);
 		Append(~M2, v);
 	end for;
 	M2:=Matrix(M2);
 	for i in [1..Nrows(M2)] do
 		for j in [1..d] do
 			M2[i][Ncols(M)+j] *:= scales[j];
 		end for;
 	end for;
 
 	ndet:=Determinant(Submatrix(M2,1,d+1,d,d));

	/* ndet is an integer because c has been chosen as (an
	 * integer)*(irrational factors that appear).
	 *
	 * If this property fails, then precision should probably be raised.
	 */
 	assert AbsoluteValue(ndet-Round(ndet)) lt 0.1;

 	if AbsoluteValue(ndet/Lbound) gt 2 then
 		printf "Uh, det=%o * Lbound\n", ndet/Lbound;
 		assert false;
 	end if;
 
 
 	M2:=Matrix(Nrows(M2),Ncols(M2),[Floor(x):x in Eltseq(M2)]);
 
 	/* The computation of the determinant seems to
 	 * encounter cancellations. The RR and CC structures have variable
 	 * precision to account for this. It is not necessarily terrible to
 	 * have things wrong here, but it's perhaps worth checking.
 	 */
 	M2:=LLL(M2);

	/* Our delta element: */
	delta:=OK!(Eltseq(M2[1])[1..d]);
	// printf "Selected %o\n", M2[1];

	if c`vals_tail eq [1,1] and #c`deltas ne 0 and s eq 1 and delta/c`deltas[#c`deltas][1] in {-1,1} then
		printf "Reached loop with delta=%o\n", delta;
		c`done:=true;
		return;
	end if;

	/*
	old:=Norm(Vector(c`lognorm));
	foo:=c`lognorm;
	foo[s] -:= Log(AbsoluteValue(Norm(delta)));
	new:=Norm(Vector(foo));

	if new gt old then
		print "Reached a point where coefficients increase";
		c`done:=true;
	end if;
	*/

	newH:=ideal<OK|delta>/I;

	if assigned(c`I) then
		assert c`I[s] subset I;
		c`I[s]/:=I;
		c`I[3-s]*:=newH;
		assert IsIntegral(c`I[s]);
		assert IsIntegral(c`I[3-s]);
	end if;

	changes:=[RR| 0, 0] ;
	/* On the canceled side, we've canceled exactly fullnorm, which
 	 * includes the old H part */
	c`H[s]:=ideal<OK|1>;
	c`lognorms[s] -:= Log(RR!fullnorm);
	changes[s] := -Log(RR!fullnorm);

	/* On the other side, we have a new H contribution */
	c`H[3-s]*:=newH;
	assert IsIntegral(newH);
	nn:=Norm(newH);

	/* This assertion is unfortunately wrong. If we do a second LLL
 	 * reduction here, then the only thing we know on the size of the
	 * vector in the end will be related to the determinant of the
	 * gram matrix of M2 ; which is a determinant of a sum of matrices.

	 * For this reason, Lbound was set above C*10^80 ; and here we check
	 * according to that.
	 */
	assert nn le C * 10^40;
	c`lognorms[3-s] +:= Log(RR!nn);
	changes[3-s] := Log(RR!nn);

	/* Concerning the logarithmic embeddings, we should have reduced the
	 * global (L2) size. Since the sum of the coefficients is constrained,
	 * this implies that the discrepancy (=unit contribution) has
	 * been reduced. */
	logemb:=sign[s] * logembeddings(delta);
	c`embs -:= logemb;
	changes cat:= Eltseq(-logemb);

	/*
	changes:=[Floor(x):x in changes];
	now:=[Floor(x):x in c`lognorms cat Eltseq(c`embs)];
	print changes, now;
	*/

	for x in selected do
		c`vals[s][x[1]]-:=x[2];
	end for;
	while c`vals_tail[s] gt 1 and c`vals[s][c`vals_tail[s]] eq 0 do
		c`vals_tail[s]-:=1;
	end while;
	Append(~c`deltas, <delta, sign[s]>);
end procedure;
