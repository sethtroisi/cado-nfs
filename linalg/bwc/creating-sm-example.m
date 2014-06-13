





if not assigned N or not assigned nsm then
    error "Must assign N and nsm first";
end if;

M:=Matrix(Integers(),N+nsm,N,[]);
for i in [1..N+nsm] do for j in [1..12] do k:=Random([1..N]); M[i,k]+:=Random(5)-2; end for; end for;

F:=Open("example.matrix.txt", "w");

fprintf F, "%o %o\n", N, N;
for i in [1..N] do
    nz:=[j:j in [1..N]|M[i,j] ne 0];
    fprintf F, "%o %o\n", #nz, Join([Sprintf("%o:%o", j-1,M[i,j]):j in nz], " ");
end for;
delete F;

F:=Open("example.matrix2.txt", "w");

fprintf F, "%o %o\n", nsm, N;
for i in [N+1..N+nsm] do
    nz:=[j:j in [1..N]|M[i,j] ne 0];
    fprintf F, "%o %o\n", #nz, Join([Sprintf("%o:%o", j-1,M[i,j]):j in nz], " ");
end for;
delete F;

sm:=[VectorSpace(GF(p),N+nsm)|];
for j in [1..nsm-1] do
    Append(~sm,Random(Universe(sm)));
end for;
/* adjust last SM so that there is a kernel mod p */
v0:=Vector(GF(p),[1..N]);

v:=v0*Transpose(Matrix(GF(p),M));
for j in [1..nsm-1] do
    v+:=(N+j)*sm[j];
end for;
Append(~sm, v/-(N+nsm));


assert Rank(Matrix(GF(p),M)) eq N;

Mfull:=HorizontalJoin(Matrix(GF(p),M),Transpose(Matrix(sm)));
assert IsZero(Vector(GF(p),[1..N+nsm])*Transpose(Mfull));

F:=Open("example.sm.txt", "w");
fprintf F, "%o\n", N;
for i in [1..N] do
    fprintf F, "%o\n", Join([IntegerToString(Integers()!sm[j][i]):j in [1..nsm]]," ");
end for;
delete F;


F:=Open("example.sm2.txt", "w");
fprintf F, "%o\n", nsm;
for i in [N+1..N+nsm] do
    fprintf F, "%o\n", Join([IntegerToString(Integers()!sm[j][i]):j in [1..nsm]]," ");
end for;
delete F;


