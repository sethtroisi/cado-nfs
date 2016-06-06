
SetSeed(1);

N:=100; nsm:=4; p:=NextPrime(2^80);
// N:=133; nsm:=5; p:=NextPrime(2^80);
// N:=133; nsm:=5; p:=1009;
N:=47; nsm:=3; p:=1009;
// N:=273; nsm:=5;  p:=NextPrime(2^130);

printf "Taken prime=%o\n", p;

function read_smat(filename)
    F:=Open(filename, "r");
    nr, nc:=Explode(StringToIntegerSequence(Gets(F)));
    M:=SparseMatrix(Integers(), nr, nc);
    for i in [0..nr-1] do
        t:=Split(Gets(F), " ");
        nt:=StringToInteger(t[1]);
        assert #t eq nt+1;
        t:=t[2..#t];
        for jj in [1..nt] do
            j,c:=Explode(Split(t[jj],":"));
            j:=StringToInteger(j);
            c:=StringToInteger(c);
            M[i+1,j+1]:=c;
        end for;
    end for;
    delete F;
    return M;
end function;

/* old */
/*
load "linalg/bwc/creating-sm-example.m";
SM_orig:=Matrix(sm);
M_orig:=M;
*/

/* new */
System(Sprintf("$PWD/build/$HOSTNAME.debug/linalg/bwc/random_matrix %o -d 12 -c 10 sm=%o,%o,$PWD/example.sm0.txt > $PWD/example.matrix0.txt", N+nsm, nsm, p));
// cut into matrix + matrix2 and sm + sm2
System(Sprintf("(echo %o %o ; tail -n +2 example.sm0.txt | head -%o) > example.sm.txt", N, nsm, N));
System(Sprintf("(echo %o %o ; tail -n +%o example.sm0.txt) > example.sm2.txt", nsm, nsm, N+2));
System(Sprintf("(echo %o %o ; tail -n +2 example.matrix0.txt | head -%o) > example.matrix.txt", N, N, N));
System(Sprintf("(echo %o %o ; tail -n +%o example.matrix0.txt) > example.matrix2.txt", nsm, N, N+2));
M_orig:=Matrix(read_smat("example.matrix0.txt"));
SM_orig:=Transpose(Matrix(GF(p),N+nsm,nsm,x[3..#x]) where x is StringToIntegerSequence(Read("example.sm0.txt")));
/*
*/





Mfull_orig:=HorizontalJoin(Matrix(GF(p),Matrix(M_orig)),Transpose(SM_orig));
assert IsZero(Vector(GF(p),[1..N+nsm])*Transpose(Mfull_orig));

MM:=read_smat("example.matrix.txt");
MM2:=read_smat("example.matrix2.txt");

SM:=Matrix(GF(p),N,nsm,x[3..#x]) where x is StringToIntegerSequence(Read("example.sm.txt"));
SM2:=Matrix(GF(p),nsm,nsm,x[3..#x]) where x is StringToIntegerSequence(Read("example.sm2.txt"));
assert SM eq Transpose(Submatrix(SM_orig,1,1,nsm,N));
assert SM2 eq Transpose(Submatrix(SM_orig,1,N+1,nsm,nsm));
assert Matrix(MM) eq Submatrix(M_orig,1,1,N,N);
Mfull:=HorizontalJoin(Matrix(GF(p),Matrix(MM)),SM);
Mfull2:=HorizontalJoin(Matrix(GF(p),Matrix(MM2)),SM2);
assert IsZero(Vector(GF(p),[1..N+nsm])*Transpose(Mfull));
assert IsZero(Vector(GF(p),[1..N+nsm])*Transpose(Mfull2));






System(Sprintf("./linalg/bwc/bwc-ptrace.sh wipe=1 prime=%o matrix=$PWD/example.matrix.txt  sm=$PWD/example.sm.txt n=6 bindir=$PWD/build/$HOSTNAME.debug/linalg/bwc", p));

load "linalg/bwc/bwc-ptrace.m";

/* Beware, at this point N has been replaced by a padded value. Recover the
 * good count.
 */
N:=Nrows(M_orig)-nsm;

/* Mx is the N*N matrix **without** SMs */






/* The following checks are only valid when Fr is the _full_ polynomial, as
 * normally computed by plingen. By taking out the rhs, we have something
 * trickier.
 */

// mmod(mdiv(A*DiagonalMatrix([X,X,X,X,1,1])*Fr,18),32);


// mmod(mdiv(A*Fr,24),76);

/*
mcoeff(A,0) eq Transpose(Matrix(n,m,&cat[Eltseq(v)[1..m]:v in RHS]));
mcoeff(A,1) eq Transpose(Matrix(n,m,&cat[Eltseq(v*Transpose(Msmall*Qsmall))[1..m]:v in RHS]));
mcoeff(A,2) eq Transpose(Matrix(n,m,&cat[Eltseq(v*Transpose(Msmall*Qsmall)^2)[1..m]:v in RHS]));
*/

// &+[mcoeff(A,i)*mcoeff(Fr,24-i):i in [0..24]];
// &+[mcoeff(A,i)*mcoeff(F,i):i in [0..24]];

/* computation organized stupidly, therefore very long */
// &+[Transpose(Matrix(4,8,&cat[Eltseq(v*Transpose(Msmall*Qsmall)^i)[1..8]:v in [y0,y1,y2,y3]]))*mcoeff(Fr,24-i):i in [0..24]];

F0:=mcoeff(F,0);
F1:=mdiv(F, 1);

/* As long as we compute at least nsm solutions out of mksol+gather, we can
 * finish the computation.  */
nsols:=#SS;
assert nsols ge #RHS;

for i in [1..nsols] do
    for j in [1..n] do
        SS[i,j] eq mpol_eval(V_0[j],Transpose(Msmall*Qsmall),F[j,i]);
    end for;
end for;

/* Evaluate all possible RHSs */

all_rhs:=Transpose(rhs)*Matrix(RHS);
ww:=[];

IsZero(all_rhs[1] + &+SS[1]*Transpose(Msmall*Qsmall));

for c in [1..nsols] do
    // v:=&+[mpol_eval(V_0[j],Transpose(Msmall*Qsmall),F[j,c]): j in [1..n]];
    v:=&+SS[c];
    IsZero(all_rhs[c] + v * Transpose(Msmall*Qsmall));

    w0:=v * Transpose(Qsmall);
    w1:=Transpose(rhs)[c];

    IsZero(w0 * Transpose(Msmall) + w1*Matrix(RHS));
    w:=Vector(Eltseq(w0)[1..N] cat Eltseq(w1));
    // IsZero(w * Transpose(Mfull));
    Append(~ww, w);
end for;

proj:=Matrix(ww) * Transpose(Mfull2);
assert Rank(proj) eq nsm-1;
comb:=Basis(Nullspace(proj))[1];

wx:=comb*Matrix(ww);
wx/:=wx[1];

print wx;

/*



sols:=Matrix(4,N,[&+[mpol_eval(Yblock[j],Transpose(Msmall*Qsmall),F[j,c]):j in [1..4]]:c in [1..4]]);


v:=Basis(Nullspace(sols))[1];

vF:=Vector(KP,v)*Transpose(F);

IsZero(&+[mpol_eval(Yblock[j],Transpose(Msmall*Qsmall),vF[j]):j in [1..4]]);

vcoeff:=func<v,j|Vector(mcoeff(Matrix(v),j))>;
vmod:=func<v,j|Vector(mmod(Matrix(v),j))>;
vdiv:=func<v,j|Vector(mdiv(Matrix(v),j))>;

vF0:=vcoeff(vF,0);
vF1:=vdiv(vF,1);
IsZero(
&+[Yblock[j]*vF0[j]:j in [1..4]]
+
&+[mpol_eval(Yblock[j],Transpose(Msmall*Qsmall),vF1[j]):j in [1..4]] * Transpose(Msmall*Qsmall));




w0:=&+[Yblock[j]*vF0[j]:j in [1..4]];
w1:=&+[mpol_eval(Yblock[j],Transpose(Msmall*Qsmall),vF1[j]):j in [1..4]];


lambda:=-vF0[4];

w1/lambda * Transpose(Q);

vF0/lambda;
*/
