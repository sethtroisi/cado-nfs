
// read the full matrix.
CLF:=Open(Sprintf("%o/complete",subdir),"r");
nerr:=0;
lv:=[RSpace(Integers(),#ideals)|];
lg:=[];

splines:=[];
while true do
	s:=Gets(CLF);
	if IsEof(s) then break; end if;
	sc:=[StringToInteger(x):x in Split(s, " :")];
	q:=sc[1];
	r:=sc[2];
	xx:=sc[3];
	assert 4+2*sc[4] eq #sc;
	v:=RSpace(Integers(),nrel)!0;

	sl:=[];

	for i in [0..sc[4]-1] do
		idx:=sc[5+2*i];
		val:=sc[6+2*i];

		v[1+idx]:=val;
		Append(~sl, <idx+1,val>);
	end for;

	v:=transform_to_ideal_valuation_vec(Eltseq(v),#ideals);
	Append(~lv,v);
	Append(~lg,<q,r,xx>);
	Append(~splines, sl);
end while;
delete CLF;

