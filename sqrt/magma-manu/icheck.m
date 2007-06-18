
/* It seems somewhat necessary to verify integrity of ideal ordering. */

print "Checking ideal ordering";
CLF:=Open(Sprintf("%o/complete",subdir),"r");
nerr:=0;
for i in [1..200] do
	s:=Gets(CLF);
	if IsEof(s) then break; end if;
	sc:=[StringToInteger(x):x in Split(s, " :")];
	q:=sc[1];
	r:=sc[2];
	xx:=sc[3];
	assert 4+2*sc[4] eq #sc;
	v:=RSpace(Integers(),nrel)!0;
	for i in [0..sc[4]-1] do
		v[1+sc[5+2*i]]:=sc[6+2*i];
	end for;
	v:=transform_to_ideal_valuation_vec(Eltseq(v),#ideals);

	w:=RSpace(Integers(),#ideals)!0;
	l:=[];
	for i in Factorization(ideal<OK|r+q*xx-a>) do
		ix:=Index(ideals,i[1]);
		w[ix]:=i[2];
		Append(~l,<ix,i[2]>);
	end for;

	if v ne w then
		l:=Sort(l, func<a,b|a[1]-b[1]>);
		me :=IntegerToString(#l);
		me cat:=&cat[Sprintf(" %o:%o",u[1],u[2]): u in l];

		ll:=[<i,v[i]>:i in [1..#ideals]|v[i] ne 0];
		prg :=IntegerToString(#ll);
		prg cat:=&cat[Sprintf(" %o:%o",u[1],u[2]): u in ll];

		nerr+:=1;
		printf "Error in ideal ordering (#%o, %o scanned):\n",nerr,i;
		print me,"\n",prg,"\n";
		if nerr gt 10 then break; end if;
	end if;
end for;
delete CLF;
if nerr ne 0 then assert false; end if;
print "ok";

