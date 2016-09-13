ps:=func<x|Join([Sprintf("%o",y):y in x], " ")>;
for i in [1..100] do
    d:=Random(10)+2;
    repeat
        f:=Polynomial([Random(-100,100):i in [0..d]]);
    until IsIrreducible(f);
    K:=NumberField(f);
    E:=Order([LeadingCoefficient(f)*K.1]);
    primes:={ pi[1]:pi in TrialDivision(Discriminant(E)) |
        pi[1] le 2^63 and Valuation(Discriminant(E),pi[1]) ge 2 };
    for p in primes do
        O:=pMaximalOrder(E,p);
        B:=BasisMatrix(O);
        Bd:=Denominator(B);
        Bz:=Matrix(Integers(),B*Bd);
        print d, ps(Eltseq(f));
        print ps([p, Bd] cat Eltseq(Bz));
    end for;
end for;
