
/* This magma script must be run from the directory containing the data files
 * from a bwc run (converted to magma format by bwc-ptrace.sh */


load "mn.m";
load "t.m";
M:=var;
M_orig:=M;

/* We'll pad the matrix to something which has nice dimensions, multiples of
 * nh*nv. Otherwise I'll get crazy very soon */
load "placemats.m";
/* Now M is M_orig zero-padded in rows & cols */
/* while Mt is the matrix as dispatched on the different threads/jobs */

KP<X>:=PolynomialRing(GF(p));

function group(nlimbs)
    return func<var|[Seqint([x mod 2^64:x in var[i*k+1..i*k+k]],2^64):i in [0..#var div k - 1]] where k is nlimbs>;
end function;

g:=group(Ceiling(Log(2^64,p)));


// assert nc le nr;
// nc:=nr;
// x:=Matrix(GF(p),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;

load "rw.m"; rw:=var;
load "cw.m"; cw:=var;

load "b.m"; 
nr_orig:=nr;
nc_orig:=nc;
nr:=tr;
nc:=tc;

/*
//L:=[]; for i in [0..nh-1] do
//t:=nr div nh + (i lt nr mod nh select 1 else 0); 
//s:=i*nrp;
//L cat:=[s+1..s+t];
//end for;Lr:=L;
//for j in [1..#Lr] do Cr[Lr[j]][j]:=1; end for;
Cr:=Matrix(GF(p),nh*nrp,nr,[]);
InsertBlock(~Cr,IdentityMatrix(GF(p), Minimum(nh*nrp,nr)),1,1);
//L:=[]; for j in [0..nv-1] do
//t:=nc  div nv + (j lt nc mod nv select 1 else 0);
//s:=j*ncp;
//L cat:=[s+1..s+t];
//end for;Lc:=L;
// for j in [1..#Lc] do Cc[Lc[j]][j]:=1; end for;
Cc:=Matrix(GF(p),nv*ncp,nc,[]);
InsertBlock(~Cc,IdentityMatrix(GF(p),Minimum(nv*ncp,nc)),1,1);

function is_almost_one(M, n)
    if not IsDiagonal(M) then return false; end if;
    for i in [1..n] do if M[i,i] ne 1 then return false; end if; end for;
    return true;
end function;

assert is_almost_one(Transpose(Cr)*Cr, Minimum(nh*nrp,nr));
assert is_almost_one(Transpose(Cc)*Cc, Minimum(nv*ncp,nc));
*/

ncx := nv*ncp;
nrx := nh*nrp;

nz:=nrp div nv;
assert nz eq ncp div nh;
pr:=func<x|(qr*nh+qq)*nz+r+1 where qq,qr is Quotrem(q, nv) where q,r is Quotrem(x-1,nz)>;
prinv:=func<x|(qr*nv+qq)*nz+r+1 where qq,qr is Quotrem(q, nh) where q,r is Quotrem(x-1,nz)>;
prp:=SymmetricGroup(nrx)![pr(x):x in [1..nrx]];
Pr:=PermutationMatrix(GF(p),[pr(x):x in [1..nrx]]);
assert Pr eq PermutationMatrix(GF(p),prp);



sc:=SymmetricGroup(nv*ncp)!colperm;
if assigned rowperm then
    sr:=SymmetricGroup(nh*nrp)!rowperm;
else
    sr:=SymmetricGroup(nh*nrp)!colperm;
end if;
sr_inv:=sr^-1;
sc_inv:=sc^-1;
Sr:=PermutationMatrix(GF(p),sr);
Sc:=PermutationMatrix(GF(p),sc);
if assigned rowperm then
    Sr:=Pr^-1*Sr;
    sr:=prp^-1*sr;
end if;

// Mx:=Cr*M*Transpose(Cc);
Mx:=M;

assert Mx[2] eq (Sr*Mx)[2^(sr^-1)];
assert Transpose(Mx)[2] eq (Sc*Transpose(Mx))[2^(sc^-1)];

S:=Sr;
P:=Pr;

// de-correlation permutation:
q:=SymmetricGroup(ncx)![preshuf(x):x in [1..ncx]];
Q:=PermutationMatrix(GF(p),q);
q0:=SymmetricGroup(nc)![preshuf(x):x in [1..nc]];
Q0:=PermutationMatrix(GF(p),q0);

/* Does this hold, really ? */
print Transpose(Pr*Sr)*Mt*Sc eq Mx*Q, " (this boolean may sometimes be wrong, I think)";
// t means twisted.

/* Note that Mt being equal to Pr times Sr*Cr*M*Transpose(Cc)*Sc^-1 is a
 * direct consequence of the fact that balancing_workhorse.c permutes the
 * rows which are being sent.
 */
if Mt eq P*S*Mx*Q*S^-1 then
    print "Shuffled product detected";
elif  Mt eq S*Mx*Q*S^-1 then
    print "non-shuffled product detected";
    Pr:=Identity(Parent(Pr));
else
    assert false;
end if;
assert Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
assert Mt eq Pr*Sr*Mx*Q*Sc^-1;

/****************************/

colweights:=func<M|[#[i:i in [1..Nrows(M)]|M[i,j] ne 0]: j in [1..Ncols(M)]]>;
rowweights:=func<M|[#[j:j in [1..Ncols(M)]|M[i,j] ne 0]: i in [1..Nrows(M)]]>;
// colsums:=func<M|[&+[Integers()!M[i,j]:i in [1..Nrows(M)]]: j in [1..Ncols(M)]]>;
// rowsums:=func<M|[&+[Integers()!M[i,j]:j in [1..Ncols(M)]]: i in [1..Nrows(M)]]>;

assert colweights(Mt*Sc) eq colweights(Mx*Q);
// assert rowweights(Transpose(Pr*Sr*Cr)*Mt) eq rowweights(M);
assert rowweights(Transpose(Pr*Sr)*Mt) eq rowweights(M);


assert IsZero((Mx*Q)[nr_orig+1..nr]);
assert IsZero(Transpose(Mx*Q)[nc_orig+1..nc]);



/* The matrix by which we multiply, really */

assert Submatrix(Transpose(Mx*Q),1,1,nr_orig,nr_orig) eq Submatrix(Transpose(Q),1,1,nr_orig,nr_orig)*Transpose(Matrix(M_orig));


Qsmall:=Submatrix(Q,1,1,nr_orig,nr_orig);
Msmall:=Matrix(GF(p),Matrix(M_orig));

MM:=Transpose(Msmall*Qsmall);


load "R.m";


// truncvec:=func<x|Vector(Eltseq(x)[1..Maximum(nr_orig, nc_orig)])>;
// expandvec:=func<x|Vector(xx cat [0:i in [#xx+1..nr]]) where xx is Eltseq(x)>;

load "VV.m";

load "C0.m"; C0:=Vector(GF(p),nr_orig,g(var));
load "C1.m"; C1:=Vector(GF(p),nr_orig,g(var));
// load "C50.m"; C50:=Vector(GF(p),nr_orig,g(var));

/* Only pick first X, since we only want to verify the check vector */
load "x.m";
Xt:=Matrix(GF(p),nr_orig,1,[]);
for i in [1..1] do for u in var[i] do Xt[u][i] +:=1; end for; end for;

Xfull:=Matrix(GF(p),nr_orig,#var,[]);
for i in [1..#var] do for u in var[i] do Xfull[u][i] +:=1; end for; end for;

x0:=Rows(Transpose(Xfull));
print "Dimension of the span of the vectors of is X", Dimension(sub<Universe(x0)|x0>);
xi:=[v*(Msmall*Qsmall) : v in x0];


// interval:=50;
interval:=1;

print "Checking consistency of check vector C1";
assert C1 eq xi[1];
// C1 eq Vector(Xt)*Transpose(Msmall*Qsmall);
print "Checking consistency of check vector C1: done";

// foo:=Vector(Xt);
// for i in [1..interval] do print i; foo *:= Transpose(Transpose(Msmall*Qsmall)); end for;
// // C50 eq foo;

print "Checking that RHS vectors are in V[1]";
for i in [1..#RHS] do
    assert RHS[i] eq VV[1][i];
end for;
print "Checking that RHS vectors are in V[1]: done";

print "Checking consistency of all vectors V_i";
for i in [1..#VV-1] do
    assert Matrix(VV[i])*Transpose(Msmall*Qsmall) eq Matrix(VV[i+1]);
end for;
print "Checking consistency of all vectors V_i: done";





function matpol_from_sequence(seq, m, n)
    assert #seq mod m*n eq 0;
    deg:=#seq div (m*n)-1;
    F:=&+[KP.1^i*      Matrix(KP,m,n,seq[i*m*n+1..(i+1)*m*n]):i in [0..deg]];
    Fr:= &+[KP.1^(deg-i)*Matrix(KP,m,n,seq[i*m*n+1..(i+1)*m*n]):i in [0..deg]];
    return F, Fr;
end function;

mcoeff:=func<M,k|
ChangeRing(Parent(M)![Coefficient(x,k):x in Eltseq(M)],BaseRing(BaseRing(M)))>;
mdiv:=func<M,k|Parent(M)![x div BaseRing(M).1^k: x in Eltseq(M)]>;
mmulx:=func<M,k|Parent(M)![x * BaseRing(M).1^k: x in Eltseq(M)]>;
mmod:=func<M,k|Parent(M)![x mod BaseRing(M).1^k: x in Eltseq(M)]>;
mval:=func<M|Minimum([Valuation(x):x in Eltseq(M)])>;
mdeg:=func<M|Maximum([Degree(x):x in Eltseq(M)])>;
mcol:=func<M,j|Transpose(M)[j]>;
mrow:=func<M,i|M[i]>;
function reciprocal(P,k)
    assert k ge Degree(P);
    return Parent(P)![Coefficient(P, k-i):i in [0..k]];
end function;

function mpol_eval(y,M,P)
    /* compute y * P(M) */
    assert CoefficientRing(P) eq CoefficientRing(M);
    s:=Parent(y)!0;
    if Degree(P) lt 0 then return s; end if;
    i:=Degree(P);
    s:=Coefficient(P, i) * y;
    while i gt 0 do
        s:=s*M;
        i:=i-1;
        s+:=Coefficient(P, i) * y;
    end while;
    return s;
end function;



load "A.m";
A:=matpol_from_sequence(g(var),m,n);

print "Checking consistency of A";
for k in [0..9] do
    assert mcoeff(A, k) eq Matrix(n,n,[(x0[i],VV[k+1][j]):i,j in [1..n]]);
end for;
print "Checking consistency of A: done";



load "F.m";
F,Fr:=matpol_from_sequence(g(var),n,n);
F:=Transpose(F);
Fr:=Transpose(Fr);
assert (1+mdeg(F))*n*n eq #var;

load "Fchunks.m";
print "Checking consistency of big file F with chunks saved in F.*.sols*";
for i,j in [1..n] do
    assert Polynomial(Fchunks[i][j]) eq F[j,i];
end for;
print "Checking consistency of big file F with chunks saved in F.*.sols*: done";

degA:=mdeg(A);
coldegs:=[mdeg(mcol(Fr,j)):j in [1..n]];

/* At this point, in the case where no "nrhs" argument has been passed to
 * plingen, the computed matrix F should be such that the following predicates
 * hold.
 */
load "rhscoeffs.m";
rhscoeffs:=Matrix(GF(p),#RHS,n,g(var));
print "Checking generator computed by plingen";

degF:=mdeg(F);

if #RHS eq 0 then
    print "(homogeneous case)";
    AF:=mmod(mdiv(A,1)*Fr,degA);
    for j in [1..n] do
        assert IsZero(mdiv(mcol(AF,j),coldegs[j]));
    end for;
else
    print "(inhomogeneous case)";
    /* In this case, we have operated on a matrix A where only the last
     * columns have been chopped off in the computation.
     */
    r:=#RHS;
    At:=Transpose(A);
    Ax:=Transpose(VerticalJoin(Matrix(At[1..r]), mdiv(Matrix(At[#RHS+1..n]),1)));
    /* And Fr is also a bit peculiar, as some rows of it have been shifted out
     * to form rhscoeffs. The rhscoeffs were in the constant coefficient in F,
     * and therefore in the head coefficients in Fr.
     */

    Ax0:=Submatrix(Ax, 1, 1, n, r);
    Ax1:=Submatrix(Ax, 1, r+1, n, n-r);
    assert Ax0 eq Submatrix(A,1,1,n,r);
    assert Ax1 eq mdiv(Submatrix(A,1,r+1, n, n-r),1);

    Fx:=VerticalJoin(rhscoeffs + mmulx(Matrix(F[1..#RHS]),1),Matrix(F[#RHS+1..n]));
    Fx0:=Submatrix(Fx, 1, 1, r, n);
    Fx1:=Submatrix(Fx, r+1, 1, n-r, n);
    assert Fx0 eq rhscoeffs + X*Submatrix(F, 1, 1, r, n);
    assert Fx1 eq Submatrix(F, r+1, 1, n-r, n);

    Frx:=Parent(Fx)![reciprocal(x,degF):x in Eltseq(Fx)];
    Frx0:=Parent(Fx0)![reciprocal(x,degF):x in Eltseq(Fx0)];
    Frx1:=Parent(Fx1)![reciprocal(x,degF):x in Eltseq(Fx1)];
    /* A has an inherent O(X^(degA+1)). Since we've shifted a few columns,
     * it's even O(X^degA).
     */
    AxFx:=mmod(Ax*Frx,degA);
    for j in [1..n] do
        assert IsZero(mdiv(mcol(AxFx,j),coldegs[j]));
    end for;

    /* Coefficient of degree coldegs[j] is column j of the product Ax*Frx */
        /* We have a sum of r contributions related to the RHS, and n-r
         * related to the random vectors */
    // Here are two ways to make the check which is made in a simpler form
    // further down. These are commented out because they really test the same
    // thing, but we keep them nevertheless because these do have some
    // interest regarding reverse engineering of the formula which is used for
    // the simpler check which comes afterwards.
    /*
    for j in [1..n] do
        d:=coldegs[j];
        assert IsZero(mcoeff(mcol(Ax0*Frx0+Ax1*Frx1,j),d));

        v:=Matrix(RHS);
        summand1:=Transpose(rhscoeffs)*Matrix(RHS);
        for k in [1..degF+1] do
            v:=v*Transpose(Msmall*Qsmall);
            assert mcoeff(Fx0,k) eq Submatrix(mcoeff(F,k-1),1,1,r, n);
            summand1+:=Transpose(mcoeff(Fx0,k))*v;
        end for;
        summand2:=Parent(summand1)!0;
        w:=Matrix(VV[1][r+1..n]);
        for k in [1..degF+1] do
            w:=w*Transpose(Msmall*Qsmall);
            assert mcoeff(Fx1,k-1) eq Submatrix(mcoeff(F,k-1),r+1, 1, n-r, n);
            summand2+:=Transpose(mcoeff(Fx1,k-1))*w;
        end for;
        IsZero(summand1 + summand2);
    end for;
    */

        /*
        summand1:=Transpose(rhscoeffs)*Matrix(RHS);
        summand2:=Parent(summand1)!0;
        // v := Matrix(VV[1][1..r]);
        // w := Matrix(VV[1][r+1..n]);
        z := Matrix(VV[1]);
        for k in [1..degF+1] do
            // v:=v*Transpose(Msmall*Qsmall);
            // w:=w*Transpose(Msmall*Qsmall);
            // assert mcoeff(Fx0,k) eq Submatrix(mcoeff(F,k-1),1,1,r, n);
            // assert mcoeff(Fx1,k-1) eq Submatrix(mcoeff(F,k-1),r+1, 1, n-r, n);
            // addend:=Transpose(mcoeff(Fx0,k))*v;
            // addend+:=Transpose(mcoeff(Fx1,k-1))*w;
            // addend eq 
            summand2+:=Transpose(mcoeff(F,k-1))*z;
            z:=z*Transpose(Msmall*Qsmall);
        end for;
        IsZero(summand1 + summand2 * Transpose(Msmall*Qsmall));
        */

end if;
print "Checking generator computed by plingen: done";


/* ok. Now think about this a bit more. This means we have many zeroes. So
 * something ought to be the null vector, right ?
 */
print "Checking that we have a solution to our system";
if #RHS gt 0 then
    images:=Transpose(rhscoeffs)*Matrix(RHS);
else
    images:=Matrix(GF(p),n, nr_orig, []);
end if;
solutions:=Parent(images)!0;
z := Matrix(VV[1]);
for k in [0..degF] do
    for j in [1..n] do
        solutions[j]+:=mcoeff(mcol(F,j),k)*z;
    end for;
    z:=z*Transpose(Msmall*Qsmall);
end for;
for j in [1..n] do
    printf "Checking that column %o (sols%o-%o) is indeed a solution... ", j, j-1, j;
    assert IsZero(images[j] + solutions[j] * Transpose(Msmall*Qsmall));
    print " ok";
end for;


load "S.m";
assert #vars mod (n*nblocks) eq 0;
nsols:=#vars div (n*nblocks);
printf "Computed %o solution vectors with mksol\n", nsols;
assert nsols ge #RHS;
SS:=[[[(Vector(GF(p),g(vars[(k*nsols+i)*n+j+1]))):j in [0..n-1]]:i in [0..nsols-1]]:k in [0..nblocks-1]];

SS:=[[&+[SS[k+1,i+1,j+1]:k in [0..nblocks-1]]:j in [0..n-1]]:i in [0..nsols-1]];

print "Checking that mksol has computed what we expect";
for i in [1..nsols] do
    for j in [1..n] do
        assert SS[i,j] eq mpol_eval(V_0[j],Transpose(Msmall*Qsmall),F[j,i]);
    end for;
end for;
print "Checking that mksol has computed what we expect: done";

all_rhs:=Transpose(rhscoeffs)*Matrix(RHS);

IsZero(all_rhs[1] + &+SS[1]*Transpose(Qsmall)*Transpose(Msmall));

// RHS has the RHS in rows.
Mpad:=HorizontalJoin(Msmall, Transpose(Matrix(RHS)));
ww:=[];
for c in [1..nsols] do
    // v:=&+[mpol_eval(V_0[j],Transpose(Msmall*Qsmall),F[j,c]): j in [1..n]];
    v:=&+SS[c];
    assert IsZero(all_rhs[c] + v * Transpose(Qsmall)*Transpose(Msmall));
    w0:=v * Transpose(Qsmall);
    w1:=Transpose(rhscoeffs)[c];
    assert IsZero(w0 * Transpose(Msmall) + w1*Matrix(RHS));
    w:=Vector(HorizontalJoin(w0,w1));
    Append(~ww, w);
end for;

print "Check that we have w such that (M||RHS) * w = 0: ",
    IsZero(Matrix(ww) * Transpose(Mpad));


print "Checking that gather has computed what we expect";
load "K.m";
assert #vars eq #RHS;
ker:=[Vector(GF(p),x):x in vars];
print "Checking that gather has computed what we expect: done";

/* FIXME -- we must get rid of Qsmall here !!!! */
print "WARNING: the data in K.sols* is not correctly de-shuffled !!";
assert IsZero(Matrix(ker)*Transpose(HorizontalJoin(Msmall*Qsmall, Transpose(Matrix(RHS)))));

exit;

// rest is utter crap.

// load "K.sols0-1.0.m";    K0:=expandvec(Vector(GF(p), g(var)));

/*
W:=K0;
for i in [1..10] do
    if IsZero(W*Transpose(Mx*Q)) then break; end if;
    printf "Now considering M^%o*W\n", i;
    W:=W*Transpose(Mx*Q);
end for;

notinker:=false;
if not IsZero(W*Transpose(Mx*Q)) then
    print "W not in the kernel of Mx*Q !?";
    System("sleep 1");
    notinker:=true;
end if;
*/
    notinker:=true;

assert Mt eq Pr*Sr*Mx*Q*Sc^-1;
print Mt eq P*S*Mx*Q*S^-1, " (true only for shuffled product)";
Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
// Mx eq Cr*M*Transpose(Cc);
Mx eq M;

// IsZero(W*Transpose(Q)*Cc*Transpose(M)*Transpose(Cr));
// sol:=Vector(Eltseq(W*Transpose(Q)*Cc)[1..nc]);













/*
if not notinker then
IsZero(W*Transpose(Q)*Transpose(M));
sol:=Vector(Eltseq(W*Transpose(Q))[1..nc]);
IsZero(sol*Transpose(M));
end if;
*/


// load "x.m";
// C0:=Matrix(GF(p),nr,64,[]);
// Xt:=Matrix(GF(p),nr,n,[]);
// for i in [1..64] do for u in var[i] do C0[u][i] +:=1; end for; end for;
// for i in [1..n] do for u in var[i] do Xt[u][i] +:=1; end for; end for;
// 
// 
// function vblock(nr,var, nc)
//     R:=Matrix(GF(p),nr,nc,[]);
//     nb:=nc div 64;
//     for i in [0..nr-1] do
//         R[i+1]:=Vector(GF(p),Intseq(u,2,nc) where u is Seqint(var[1+i*nb..(i+1)*nb],2^64));
//     end for;
//   return R;
// end function;
// 
// 
// load "Y.0.m";      Y0:= vblock(nc, var, n);
// load "V0-64.0.m";  V0:= vblock(nc, var, 64);
// // load "V0-64.1.m";  V1:= vblock(nc, var, 64);
// // load "V0-64.2.m";  V2:= vblock(nc, var, 64);
// // load "V0-64.3.m";  V3:= vblock(nc, var, 64);
// // load "V0-64.4.m";  V4:= vblock(nc, var, 64);
// // load "V0-64.5.m";  V5:= vblock(nc, var, 64);
// // load "V0-64.6.m";  V6:= vblock(nc, var, 64);
// // load "V0-64.7.m";  V7:= vblock(nc, var, 64);
// // load "V0-64.8.m";  V8:= vblock(nc, var, 64);
// // load "V0-64.9.m";  V9:= vblock(nc, var, 64);
// load "V0-64.10.m"; V10:=vblock(nc, var, 64);
// load "V0-64.20.m"; V20:=vblock(nc, var, 64);
// 
// load "V64-128.0.m";  V0b:= vblock(nc, var, 64);
// load "V64-128.10.m"; V10b:=vblock(nc, var, 64);
// load "V64-128.20.m"; V20b:=vblock(nc, var, 64);
// 
// 
// // load "V0-64.12.m"; V12:=vblock(nc, var, n);
// // load "V0-64.13.m"; V13:=vblock(nc, var, n);
// // load "V0-64.14.m"; V14:=vblock(nc, var, n);
// // load "V0-64.15.m"; V15:=vblock(nc, var, n);
// // load "V0-64.16.m"; V16:=vblock(nc, var, n);
// // load "V0-64.17.m"; V17:=vblock(nc, var, n);
// // load "V0-64.18.m"; V18:=vblock(nc, var, n);
// load "C.0.m";    C0:=vblock(nr, var, 64);
// load "C.10.m";   C10:=vblock(nr, var, 64);
// 
// load "H1.m";    H1:=vblock(nc, var, 64);
// var:=[(p div j+q*j) mod 2^64:j in [1..nc]]
//     where p is PreviousPrime(2^64)
//     where q is PreviousPrime(2^63);
// H0:=vblock(nc, var, 64);
// assert M*Q0*H0 eq H1;
// 
// // v0z:=[Seqint(ChangeUniverse((Eltseq(V0[i])),Integers()),2):i in [1..nr]];
// // v0bz:=[Seqint(ChangeUniverse((Eltseq(V0b[i])),Integers()),2):i in [1..nr]];
// 
// // it's somewhat misleading sometimes.
// // [Index(Rows(Mt*V0),r):r in Rows(V1)];
// 
// Z:=Integers();
// 
// /*
// Sort(rowweights(Matrix(Z,Transpose(Mt)*V0))) eq 
// Sort(rowweights(Matrix(Z,V1)));
// 
// ([w[pr(i)]: i in [1..#w]] where w is rowweights(Matrix(Z,Transpose(Mt)*V0)))\
//  eq rowweights(Matrix(Z,V1));
// 
// rowweights(Matrix(Z,Transpose(Mt)*V0)) eq
// [w[prinv(i)]: i in [1..#w]] where w is rowweights(Matrix(Z,V1));
// 
// rowweights(Matrix(Z,Pr*Transpose(Mt)*V0)) eq rowweights(Matrix(Z,V1));;
// */
// 
// TMt:=xtr(Pr^-1)*xtr(Mt);
// 
// TM:=xtr(M*Q0);
// 
// // equivalent to TMt:=Transpose(Mt*Pr^-1) for nullspace=left
// // For nullspace=left:
// // TMt eq Transpose(Mt*Pr^-1);
// // TMt eq Transpose(Pr*Sr*Cr*M*Transpose(Cc)*Sc^-1*Pr^-1);
// // TMt eq Transpose(Pr*Sr*Mx*Sc^-1*Pr^-1);
// // TMt eq Transpose(Pr*Sc)^-1*Transpose(Mx)*Transpose(Pr*Sr);
// // assert Sr eq Sc
// // Transpose(Mx)*(Pr*Sr)^-1*V0 eq (Pr*Sr)^-1*V1;
// //
// // For nullspace=right
// // TMt eq Pr^-1*Mt
// // TMt eq Sr*Cr*M*Transpose(Cc)*Sc^-1
// // TMt eq Sr*Mx*Sc^-1
// // assert Sr eq Sc
// // Mx*Sr^-1*V0 eq Sr^-1*V1
// 
// conj_if_left:=func<x|xtr(xtr(x)*Pr^-1)*xtr(Pr)>;
// // for nullspace == left this is Pr*x*Pr^-1
// 
// // it's normal if one of the two is false. It would help chasing bugs, in
// // fact !
// TMt eq conj_if_left(Sc*Cc*TM*Transpose(Cr)*Sr^-1);
// 
// // Mtt eq Pr*Sc*Cr*M*Transpose(Cc)*Sr^-1*Pr^-1;
// 
// TM*V0 eq V1; assert TM*V0 eq V1;
// TM*V1 eq V2; assert TM*V1 eq V2;
// TM*V2 eq V3; assert TM*V2 eq V3;
// TM*V3 eq V4; assert TM*V3 eq V4;
// TM^10*V0 eq V10; assert TM^10*V0 eq V10;
// TM^10*V0b eq V10b; assert TM^10*V0b eq V10b;
// 
// // (Pr^-1*Mt)^10*C0 eq C10;
// //
// // This one only valid for nullspace==left
// // Pr*(Pr^-1*Mt)^10*Pr^-1*C0 eq C10;
// //
// print "check vectors:";
// C0 eq Submatrix(Xt,1,1,Nrows(Xt),64); 
// C10 eq Submatrix(Transpose(TM)^10*Xt,1,1,Nrows(Xt),64);
// 
// Transpose(TM)^10*C0 eq C10;
// 
// 
// Transpose(V0)*C10 eq Transpose(V10)*C0;
// Transpose(C10)*V0 eq Transpose(C0)*V10;
// Transpose(C10)*V0b eq Transpose(C0)*V10b;
// 
// print "Checking krylov stuff: done";
// 
// load "A0-128.0-30.m";   // suited for matrix t100b
// 
// stride:=m*n div 64;
// for i in [0..30-1] do
//     vblock(m,var[i*stride+1..i*stride+stride],n) eq
//     Transpose(Xt)*TM^i*Y0;
// end for;
// 
K:=GF(p);
b:=m+n;
N:=nr;
// KR<x>:=PowerSeriesRing(K);
KR<x>:=PolynomialAlgebra(K);
KRmn:=RMatrixSpace(KR,m,n);
KRbb:=RMatrixSpace(KR,b,b);
KRmb:=RMatrixSpace(KR,m,b);
KRnb:=RMatrixSpace(KR,n,b);
KRn1:=RMatrixSpace(KR,n,1);
KRnn:=RMatrixSpace(KR,n,n);
KRmn:=RMatrixSpace(KR,m,n);
KRNm:=RMatrixSpace(KR,N,m);
KRNn:=RMatrixSpace(KR,N,n);
KRN:=RMatrixSpace(KR,N,N);
KNn:=KMatrixSpace(K,N,n);


// 
// // beware -- we're using the sequence starting at 1.
// // A_sequence:=&+[x^(i-1)*KRmn!(Transpose(Xt)*TM^i*V0):i in [1..100]];
// A_sequence:=KRmn!0;
// z:=Y0;
// for i in [1..100] do
//     A_sequence +:= x^(i-1)*KRmn!(Transpose(Xt)*z);
//     z:=TM*z;
// end for;
// // assert A_sequence eq &+[x^(i-1)*KRmn!(Transpose(Xt)*TM^(i-1)*V0):i in [1..100]];
// 
// 
// load "F.m";
// 
// degf:=#var div stride - 1;
// 
// F:=&+[x^(degf-i)*KRnn!vblock(n,var[i*stride+1..(i+1)*stride],n):i in [0..degf]];
// Fr:=&+[x^i*      KRnn!vblock(n,var[i*stride+1..(i+1)*stride],n):i in [0..degf]];
// 
// printf "FYE: %o\n", IsZero(mycoeff(A_sequence*Transpose(F),degf));
// for i in [degf+1..degf+5] cat [97..99] do
//     // degf or degf+1 ? I'm seeing non-null at degf... It's pretty harmless.
//     IsZero(mycoeff(A_sequence*Transpose(F),i));
// end for;
// 
// not IsZero(mycoeff(A_sequence*Transpose(F),degf-1));
// 
// excess:=10;
// 
// V_sequence:=&+[x^i*KRNn!(TM^i*Y0):i in [0..degf+excess]];
// 
// assert IsZero(Transpose(Xt)*mycoeff(V_sequence*Transpose(F),degf+1));
// 
// // Transpose(Xt)*(KRN!TM)^4*(KRN!TM*V_sequence) eq mycoeff(A_sequence,3);
// 
// IsZero(mycoeff(Transpose(KRNm!Xt)*(&+[x^(i-1)*KRNn!(TM^(i-1)*Y0):i in [1..10\
// 0]])*Transpose(F),degf+1));
// 
// 
// for i in [0..excess-1] do
//     IsZero(Transpose(Xt)*TM*mycoeff(V_sequence*Transpose(F),degf + i));
//     IsZero(TM*mycoeff(V_sequence*Transpose(F),degf + i));
// end for;
// 
// S:=&+[TM^i*Y0*mycoeff(Transpose(Fr),i):i in [0..degf]];
// IsZero(TM*S);
// 
// // S:=&+[TM^i*V0*mycoeff(Transpose(F),i):i in [0..2]];
// 
// load "S0-64.10.m"; S10:=vblock(nc, var, n);
// load "S0-64.20.m"; S20:=vblock(nc, var, n);
// load "S64-128.10.m"; S10b:=vblock(nc, var, n);
// load "S64-128.20.m"; S20b:=vblock(nc, var, n);
// load "K.0.m";    K0:=vblock(nc, var, 64);
// 
// S eq S10+S20;
// 
// // K0 is basically S in column echelon form.
// Rowspace(Transpose(S)) eq Rowspace(Transpose(K0));
// IsZero(TM*K0);
// 
// 
// 
// load "W.m";    W:=vblock(nr, var, 64);    
// 
// 
// // For our simple test cases, K0==W. In more subtle situations, things
// // might differ a bit.
// 
// print "FYI: ", K0 eq W;
// IsZero(TM*W);
