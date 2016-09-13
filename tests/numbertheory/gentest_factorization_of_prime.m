ps:=func<x|Join([Sprintf("%o",y):y in x], " ")>;
for i in [1..100] do
    d:=Random(10)+2;
    repeat
        f:=Polynomial([Random(-100,100):i in [0..d]]);
    until IsIrreducible(f) and Degree(f) ge 2;
    disc:=LeadingCoefficient(f)*Discriminant(f);
    primes:={ pi[1]:pi in TrialDivision(disc) |
        pi[1] le 2^63 and Valuation(disc,pi[1]) ge 2 };
    primes join:={2,3,5,7};
    for p in primes do
        K:=NumberField(f);
        E:=Order([LeadingCoefficient(f)*K.1]);
        O:=pMaximalOrder(E, p);
        B:=BasisMatrix(O);
        Bd:=Denominator(B);
        Bz:=Matrix(Integers(),B*Bd);
        O`Maximal:=true;
        Fp:=Factorization(ideal<O|p>);
        print Degree(f), ps(Eltseq(f));
        iprint:=function(Im)
            I,m:=Explode(Im);
            MI:=BasisMatrix(I);
            return Seqlist(Eltseq(BasisMatrix(I))) cat [* "", m *];
        end function;
        print ps([p]),
            ps([Bd] cat Eltseq(Bz)),
            ps([""]),
            ps([*#Fp*] cat &cat [[*""*] cat iprint(Im):Im in Fp]);
    end for;
end for;
