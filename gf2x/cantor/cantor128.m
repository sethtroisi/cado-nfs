/* Magma code that checks the cantor_main output: */

// w := 64;

load "/tmp/toto";
PP := PolynomialRing(GF(2));
F := PP!Intseq(Seqint(f, 2^w), 2);
G := PP!Intseq(Seqint(g, 2^w), 2);
FG := PP!Intseq(Seqint(fg, 2^w), 2);

FG eq F*G;
if FG ne F*G then
    exit;
end if;


// w := 64;
load "/tmp/toto";
PP := PolynomialRing(GF(2));
x := PP.1;
F128<z> := ext<GF(2) | x^128 + x^7 + x^2 + x + 1>;
P128<x> := PolynomialRing(F128);


function convertF128(F128, x)
  return F128!PP!Intseq(x, 2);
end function;

function convertP128(F128, f)
  return Polynomial([ convertF128(F128, z): z in Intseq(Seqint(f, 2^w), 2^64)]);
end function;

ffi := convertP128(P128,F128,fi);
ggi := convertP128(P128,F128,gi);
h := convertP128(P128,F128,fgi);
fg := ffi*ggi;


beta0 := F128!1;
beta1 := F128!Intseq((7451958335100816136 +2^64*2979905974001933049 ), 2);
beta2 := F128!Intseq((18379359142562139938+2^64*11678753174964950217), 2);
beta3 := F128!Intseq((6032672185962850376 +2^64*8744256601846606146 ), 2);
beta4 := F128!Intseq((14608789766813931076+2^64*645900494672607458  ), 2);
beta5 := F128!Intseq((5568895992531017964 +2^64*5316835906050072984 ), 2);
beta6 := F128!Intseq((11619261390503595532+2^64*377604546732988956  ), 2);
beta7 := F128!Intseq((6679137017075335448 +2^64*5571574281931689094 ), 2);

s1:=x^2+x;
s2:=Evaluate(s1,s1);
s3:=Evaluate(s2,s1);
s4:=Evaluate(s2,s2);
s5:=Evaluate(s1,s4);
s6:=Evaluate(s3,s3);
s7:=Evaluate(s3,s4);


function split(fi)
    return [Intseq(Seqint(ChangeUniverse(Eltseq(Coefficient(fi,k)),Integers()),2),2^w) : k in [0..Degree(fi)] ];
end function;

///////////////////////////////////////////////////////////////////////////

procedure quotremSi(i, pf, ii)
  printf "QuotremSi(%o, %o, %o)\n", i, pf, ii;
end procedure;

procedure mulSi(i, pf, ii)
  printf "mulSi(%o, %o, %o)\n", i, pf, ii;
end procedure;



procedure reduceModTrunc(k, length)
  pf := 0;
  len := length;

  // go down, computing quotrems
  printf "Going down...\n";
  ii := 2^k;
  for i := k-1 to 0 by -1 do
    if (len ge (2^i)) then
      quotremSi(i, pf, ii);
      ii -:= 2^i;
      pf +:= 2^i;
      len -:= 2^i;
    end if;
  end for;

  printf "Going up...\n";
  len := length;
  i := 0;
  while IsEven(len) do
    len div:= 2;
    i +:= 1;
  end while;
  pf := length;
  ii := 2^i;

  i +:= 1;
  len div:= 2;
  while i le k-1 do
    if ( (len mod 2) eq 1) then
      mulSi(i, pf - ii, ii);
      ii +:= 2^i;
    end if;
    len div:= 2;
    i +:= 1;
  end while;
end procedure;

