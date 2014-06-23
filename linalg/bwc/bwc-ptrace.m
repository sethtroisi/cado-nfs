


load "/tmp/bwcp/mn.m";
load "/tmp/bwcp/t.m";
M:=var;

/* We'll pad the matrix to something which has nice dimensions, multiples of
 * nh*nv. Otherwise I'll get crazy very soon
 */
load "/tmp/bwcp/placemats.m";

KP<X>:=PolynomialRing(GF(p));


// assert nc le nr;
// nc:=nr;
// x:=Matrix(GF(p),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;

load "/tmp/bwcp/rw.m"; rw:=var;
load "/tmp/bwcp/cw.m"; cw:=var;

load "/tmp/bwcp/b.m"; 
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

Mx[2] eq (Sr*Mx)[2^(sr^-1)];
Transpose(Mx)[2] eq (Sc*Transpose(Mx))[2^(sc^-1)];

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



function group(nlimbs)
    return func<var|[Seqint([x mod 2^64:x in var[i*k+1..i*k+k]],2^64):i in [0..#var div k - 1]] where k is nlimbs>;
end function;

g:=group(Ceiling(Log(2^64,p)));

/* No load in conditionals. Anyway we never use Y.
if n eq 1 then
load "/tmp/bwcp/Y0.m"; Y0:=Vector(GF(p),nr,g(var));
end if;
*/

truncvec:=func<x|Vector(Eltseq(x)[1..Maximum(nr_orig, nc_orig)])>;
expandvec:=func<x|Vector(xx cat [0:i in [#xx+1..nr]]) where xx is Eltseq(x)>;

load "/tmp/bwcp/V0.0.m"; V0:=expandvec(Vector(GF(p),g(var)));
load "/tmp/bwcp/V0.1.m"; V1:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.2.m"; V2:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.3.m"; V3:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.4.m"; V4:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.5.m"; V5:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.6.m"; V6:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.7.m"; V7:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.8.m"; V8:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.9.m"; V9:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/V0.10.m"; V10:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/C0.m"; C0:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/C1.m"; C1:=expandvec(Vector(GF(p),nr_orig,g(var)));
load "/tmp/bwcp/C50.m"; C50:=expandvec(Vector(GF(p),nr_orig,g(var)));

/* Only pick first X, since we only want to verify the check vector */
load "/tmp/bwcp/x.m";
Xt:=Matrix(GF(p),nr,1,[]);
for i in [1..1] do for u in var[i] do Xt[u][i] +:=1; end for; end for;


interval:=50;

C1 eq Vector(Xt)*Mx*Q;

foo:=Vector(Xt);
for i in [1..50] do print i; foo *:= Mx*Q; end for;
C50 eq foo;


V1  eq V0*Transpose(Mx*Q);
V2  eq V1*Transpose(Mx*Q);
V3  eq V2*Transpose(Mx*Q);
V4  eq V3*Transpose(Mx*Q);
V5  eq V4*Transpose(Mx*Q);
V6  eq V5*Transpose(Mx*Q);
V7  eq V6*Transpose(Mx*Q);
V8  eq V7*Transpose(Mx*Q);
V9  eq V8*Transpose(Mx*Q);
V10 eq V9*Transpose(Mx*Q);

function matpol_from_sequence(seq, m, n)
    assert #seq mod m*n eq 0;
    deg:=#seq div (m*n)-1;
    F:=&+[KP.1^i*      Matrix(KP,m,n,seq[i*m*n+1..(i+1)*m*n]):i in [0..deg]];
    Fr:= &+[KP.1^(deg-i)*Matrix(KP,m,n,seq[i*m*n+1..(i+1)*m*n]):i in [0..deg]];
    return F, Fr;
end function;

mcoeff:=func<M,k|Matrix(CoefficientRing(CoefficientRing(M)),Nrows(M),Ncols(M),[Coefficient(M[i,j],k):j in [1..Ncols(M)], i in [1..Nrows(M)]])>;
mdiv:=func<M,k|Matrix(KP,Nrows(M),Ncols(M),[M[i,j] div KP.1^k:j in [1..Ncols(M)], i in [1..Nrows(M)]] where KP is CoefficientRing(M))>;
mmod:=func<M,k|Matrix(KP,Nrows(M),Ncols(M),[M[i,j] mod KP.1^k:j in [1..Ncols(M)], i in [1..Nrows(M)]] where KP is CoefficientRing(M))>;
mval:=func<M|Minimum([Valuation(M[i,j]):i in [1..Nrows(M)], j in [1..Ncols(M)]])>;


load "/tmp/bwcp/A.m";
A:=matpol_from_sequence(g(var),m,n);


load "/tmp/bwcp/F.m";
F,Fr:=matpol_from_sequence(g(var),n,n);
F:=Transpose(F);
Fr:=Transpose(Fr);

load "/tmp/bwcp/rhs.m";
rhs:=Matrix(GF(p),nsm,n,g(var));

degF:=Maximum([Degree(F[i,j]):i,j in [1..n]]);

if IsZero(mmod(F, 1)) then
    print "Very weird, first coefficient in F stored as zero !?!!?";
end if;

for j in [1..n] do
    for e in [0..10] do
        degcol:=Maximum([Degree(Fr[i,j]):i in [1..n]]);
        v:=mval(mdiv(A*Submatrix(Fr,1,j,n,1),degcol+e));
        if v ge 10 then
            printf "Solution column %o: with delta=%o+%o, found %o zeros\n",
                j, degcol, e, v;
            break;
        end if;
    end for;
end for;
// something like this should work.
// IsZero(mmod(mdiv(A*Fr,degF),159));

load "/tmp/bwcp/V0.0.m"; y0:=expandvec(Vector(GF(p),g(var)));
load "/tmp/bwcp/V1.0.m"; y1:=expandvec(Vector(GF(p),g(var)));
load "/tmp/bwcp/V2.0.m"; y2:=expandvec(Vector(GF(p),g(var)));
load "/tmp/bwcp/V3.0.m"; y3:=expandvec(Vector(GF(p),g(var)));
load "/tmp/bwcp/V4.0.m"; y4:=expandvec(Vector(GF(p),g(var)));
load "/tmp/bwcp/V5.0.m"; y5:=expandvec(Vector(GF(p),g(var)));


load "/tmp/bwcp/S.m";
nsols:=#vars div n;
printf "Computed %o solution vectors with mksol\n", nsols;
assert #vars mod n eq 0;
SS:=[[expandvec(Vector(GF(p),g(vars[i*n+j+1]))):j in [0..n-1]]:i in [0..nsols-1]];

load "/tmp/bwcp/K.sols0-1.0.m";    K0:=expandvec(Vector(GF(p), g(var)));

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

assert Mt eq Pr*Sr*Mx*Q*Sc^-1;
print Mt eq P*S*Mx*Q*S^-1, " (true only for shuffled product)";
Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
// Mx eq Cr*M*Transpose(Cc);
Mx eq M;

// IsZero(W*Transpose(Q)*Cc*Transpose(M)*Transpose(Cr));
// sol:=Vector(Eltseq(W*Transpose(Q)*Cc)[1..nc]);
if not notinker then
IsZero(W*Transpose(Q)*Transpose(M));
sol:=Vector(Eltseq(W*Transpose(Q))[1..nc]);
IsZero(sol*Transpose(M));
end if;


// load "/tmp/bwcp/x.m";
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
// load "/tmp/bwcp/Y.0.m";      Y0:= vblock(nc, var, n);
// load "/tmp/bwcp/V0-64.0.m";  V0:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.1.m";  V1:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.2.m";  V2:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.3.m";  V3:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.4.m";  V4:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.5.m";  V5:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.6.m";  V6:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.7.m";  V7:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.8.m";  V8:= vblock(nc, var, 64);
// // load "/tmp/bwcp/V0-64.9.m";  V9:= vblock(nc, var, 64);
// load "/tmp/bwcp/V0-64.10.m"; V10:=vblock(nc, var, 64);
// load "/tmp/bwcp/V0-64.20.m"; V20:=vblock(nc, var, 64);
// 
// load "/tmp/bwcp/V64-128.0.m";  V0b:= vblock(nc, var, 64);
// load "/tmp/bwcp/V64-128.10.m"; V10b:=vblock(nc, var, 64);
// load "/tmp/bwcp/V64-128.20.m"; V20b:=vblock(nc, var, 64);
// 
// 
// // load "/tmp/bwcp/V0-64.12.m"; V12:=vblock(nc, var, n);
// // load "/tmp/bwcp/V0-64.13.m"; V13:=vblock(nc, var, n);
// // load "/tmp/bwcp/V0-64.14.m"; V14:=vblock(nc, var, n);
// // load "/tmp/bwcp/V0-64.15.m"; V15:=vblock(nc, var, n);
// // load "/tmp/bwcp/V0-64.16.m"; V16:=vblock(nc, var, n);
// // load "/tmp/bwcp/V0-64.17.m"; V17:=vblock(nc, var, n);
// // load "/tmp/bwcp/V0-64.18.m"; V18:=vblock(nc, var, n);
// load "/tmp/bwcp/C.0.m";    C0:=vblock(nr, var, 64);
// load "/tmp/bwcp/C.10.m";   C10:=vblock(nr, var, 64);
// 
// load "/tmp/bwcp/H1.m";    H1:=vblock(nc, var, 64);
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
// print "Done checking krylov stuff";
// 
// load "/tmp/bwcp/A0-128.0-30.m";   // suited for matrix t100b
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

function mycoeff(M,k)
        s:=KMatrixSpace(CoefficientRing(CoefficientRing(M)),Nrows(M), Ncols(M));
        return s![Coefficient(P,k):P in Eltseq(M)];
end function;
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
// load "/tmp/bwcp/F.m";
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
// load "/tmp/bwcp/S0-64.10.m"; S10:=vblock(nc, var, n);
// load "/tmp/bwcp/S0-64.20.m"; S20:=vblock(nc, var, n);
// load "/tmp/bwcp/S64-128.10.m"; S10b:=vblock(nc, var, n);
// load "/tmp/bwcp/S64-128.20.m"; S20b:=vblock(nc, var, n);
// load "/tmp/bwcp/K.0.m";    K0:=vblock(nc, var, 64);
// 
// S eq S10+S20;
// 
// // K0 is basically S in column echelon form.
// Rowspace(Transpose(S)) eq Rowspace(Transpose(K0));
// IsZero(TM*K0);
// 
// 
// 
// load "/tmp/bwcp/W.m";    W:=vblock(nr, var, 64);    
// 
// 
// // For our simple test cases, K0==W. In more subtle situations, things
// // might differ a bit.
// 
// print "FYI: ", K0 eq W;
// IsZero(TM*W);
