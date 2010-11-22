// Magma code called by check_cantor.sh to check the result computed
// by cantor_main.
//

load "/tmp/toto";
PP := PolynomialRing(GF(2));
x := PP.1;
F128<z> := ext<GF(2) | x^128 + x^7 + x^2 + x + 1>;
P128<x> := PolynomialRing(F128);
F := PP!Intseq(Seqint(f, 2^w), 2);
G := PP!Intseq(Seqint(g, 2^w), 2);
FG := PP!Intseq(Seqint(fg, 2^w), 2);

FG eq F*G;
