/*
vals:=[<i,vals[i]> : i in [1..nrel] | vals[i] ne 0];
vals_plus:=[v:v in vals|v[2] gt 0];
vals_minus:=[<v[1],-v[2]>:v in vals|v[2] lt 0];
vals:=[vals_plus, vals_minus];
*/

/* Hplus is H[1], Hminus is H[2] */

lengthf:=Sqrt(Norm(Vector(Eltseq(pol))));
C:=lengthf^(d-1)*2^(d*(d-1)/4);
Lbound:=Maximum(C^2*10^80,10^400);

sq_desc:=recformat<
	deltas : SeqEnum,
	H :  SeqEnum,
	I : SeqEnum,
	vals : SeqEnum,
	vals_tail : SeqEnum,
	lognorms : SeqEnum,
	embs : ModTupFldElt,
	done : BoolElt >;

c:=rec<sq_desc | 
	deltas:=[],
	H:=[ideal<OK|1>, ideal<OK|1>],
	vals:=vals0,
	vals_tail:=[#ideals,#ideals],
	lognorms := lognorms0,
	embs:=embs0 >;

c0:=c;

	/*
c`lvecs[1][1]:=Floor(&+[Log(Norm(ideals[i]))*vals0[1][i]:i in [1..#ideals]]);
c`lvecs[2][1]:=Floor(&+[Log(Norm(ideals[i]))*vals0[2][i]:i in [1..#ideals]]);
c`lvecs[1][3]:=c`lvecs[1][1];
c`lvecs[2][3]:=c`lvecs[2][1];
*/

function pick_ideals(s,c)
	selected:={};
	fullnorm:=Norm(c`H[s]);
	pos:=c`vals_tail[s];
	while fullnorm lt Lbound do
		if pos eq 0 then break; end if;
		if c`vals[s][pos] eq 0 then pos-:=1; continue; end if;
		ni:=Norm(ideals[pos]);
		maxi:=Floor(Log(Lbound/fullnorm)/Log(ni));
		maxi:=Minimum([c`vals[s][pos], maxi]);
		Include(~selected, <pos,maxi>);
		fullnorm*:=Norm(ideals[pos])^maxi;
		assert fullnorm lt Lbound;
		if maxi ne c`vals[s][pos] then break; end if;
		pos-:=1;
	end while;
	return selected, fullnorm;
end function;


	/*
procedure recompute_logs(~c)
	newlogs:=[
		[c`lvecs[1][1],
		&+[RR|Log(AbsoluteValue(Norm(d[1])))*d[2] : d in c`deltas | d[2] gt 0],
		&+[Log(Norm(ideals[i]))*c`vals[1][i]:i in [1..#ideals]],
		Log(Norm(c`H[1]))],
		[c`lvecs[2][1],
		&+[RR|-Log(AbsoluteValue(Norm(d[1])))*d[2] : d in c`deltas | d[2] lt 0],
		&+[Log(Norm(ideals[i]))*c`vals[2][i]:i in [1..#ideals]],
		Log(Norm(c`H[2]))]];

	newlogs:=[[Floor(y):y in x]:x in newlogs];
	dif:=[[newlogs[j][i]-c`lvecs[j][i] : i in [1..4]]:j in [1..2]];
	c`lvecs:=newlogs;
	// print Matrix(d);
end procedure;

function check(c)
	lhs_plus:= &+[Log(Norm(ideals[i]))*vals0[1][i] :i in [1..#ideals]];
	lhs_minus:=&+[Log(Norm(ideals[i]))*vals0[2][i] :i in [1..#ideals]];
	rhs_plus:= &+[Log(Norm(ideals[i]))*c`vals[1][i]:i in [1..#ideals]];
	rhs_minus:=&+[Log(Norm(ideals[i]))*c`vals[2][i]:i in [1..#ideals]];
	rhs_plus+:= &+[RR|Log(AbsoluteValue(Norm(d[1])))*d[2] : d in c`deltas | d[2] gt 0];
	rhs_minus+:=&+[RR|-Log(AbsoluteValue(Norm(d[1])))*d[2] : d in c`deltas | d[2] lt 0];
	rhs_plus+:= Log(Norm(c`H[1]));
	rhs_minus+:=Log(Norm(c`H[2]));

	return (lhs_plus-lhs_minus)-(rhs_plus - rhs_minus);
end function;
*/



