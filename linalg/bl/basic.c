#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "alloc.h"
#include "random.h"
#include "mat_ops.h"
#include "echelon.h"
#include "basic.h"
#include "stdio.h"

#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))


////////////////////////////////////////////////////////
//
//
//  Functions used in Lanczos 
//
////////////////////////////////////////////////////////



unsigned long Refresh_ArrayBit(unsigned long m, unsigned long n,
			       unsigned long **A, unsigned long *a)
{
    unsigned long j;
    // if (n%N==0) SizeA=n/N*m; else SizeA=n/N*m+m;
    for (j = 0; j < iceildiv(n, WBITS) * m; ++j) {
	A[2][j] = A[1][j];
	A[1][j] = A[0][j];
	A[0][j] = a[j];
    }
    return 0;
}






unsigned long TestZero(unsigned long m, unsigned long n, unsigned long *a)
{
    unsigned long i;
    for (i = 0; i < iceildiv(n, WBITS) * m; ++i)
	if (a[i] != 0)
	    return 0;
    return 1;

}




void InnerProducts(unsigned long m, unsigned long n, unsigned long N,
		   unsigned long *a, unsigned long *v, unsigned long *AV,
		   unsigned long *VAV, unsigned long *VA2V)
{
    unsigned long *ResmN;
    ResmN = Allocmn(m, N);

    // SMultDmatrixBit(m,n,N,a,v,ResmN);
    // TSMultDmatrixBit(m,n,N,a,ResmN,AV); 

    unsigned long *ResNn;
    ResNn = Allocmn(N, n);

    TransposeBit(n, N, v, ResNn);
    DMultBit(N, n, N, ResNn, AV, VAV);
    TransposeBit(n, N, AV, ResNn);
    DMultBit(N, n, N, ResNn, AV, VA2V);


    free(ResNn);
    free(ResmN);
}





char * TimeConvert(float a)
{ 

char *f;
f=malloc(10*sizeof(unsigned long));

unsigned long years,minutes,hours,days;
float seconds;

if (a>=60) {  

	if ((a/60)>=60) {

        	if ((a/3600)>=24) {
			if (a/(86400)>=365) {
					years=a/31536000;
                                        days=a/86400-years*365;
                        		hours=a/3600-years*365*24-days*24;
                        		minutes=a/60-years*365*24*60-days*24*60-hours*60;
					seconds=a-years*365*86400-days*86400-hours*3600-minutes*60;
					sprintf(f, "%luy%lud%luh%lum%.3fs",years,days,hours,minutes,seconds);
				
			}
			else
			{
                        days=a/(86400);
                        hours=a/3600-days*24;
                        minutes=a/60-days*24*60-hours*60;
			seconds=a-days*86400-hours*3600-minutes*60;
			sprintf(f, "%lud%luh%lum%.3fs",days,hours,minutes,seconds);
			}		


		}
		else
		{
                 hours=a/3600;
                 minutes=a/60-hours*60;
		 seconds=a-hours*3600-minutes*60;
		 sprintf(f, "%luh%lum%.3fs",hours,minutes,seconds);
		}

	}
	else
	{
        minutes=a/60;
        seconds=a-minutes*60;
	sprintf(f, "%lum%.3fs",minutes,seconds);
	}



}

else

{

sprintf(f, "%.3fs",a);

 
}

return f;

}







////////////////////////////////////////////////////
//
//
//  The Lanczos Algorithm  (Faster Version with a level two recursion)
//  With good multiplication.
//  Final Version
//
/////////////////////////////////////////////////////



unsigned long LanczosOld(unsigned long m, unsigned long n, unsigned long N,
			 unsigned long *a, unsigned long *v, unsigned long *X)
{
    unsigned long i, SizeRowV;

    //initial conditions data prodution

    SizeRowV = iceildiv(N, WBITS);

    unsigned long *IN;

    IN = Allocmn(N, N);

    memset(IN, 0, N * SizeRowV * sizeof(unsigned long));

    for (i = 0; i < N; ++i) {
	IN[i * SizeRowV + i / WBITS] = (1UL << i % WBITS);
    }


    unsigned long *ZN;
    ZN = Allocmn(N, N);

    memset(ZN, 0, N * SizeRowV * sizeof(unsigned long));

    unsigned long *ZnN;
    ZnN = Allocmn(n, N);

    memset(ZnN, 0, n * SizeRowV * sizeof(unsigned long));

    unsigned long **V;
    V = Allocmn3(n, N);



    for (i = 0; i < n * SizeRowV; ++i) {
	V[0][i] = v[i];
	V[1][i] = ZnN[i];
	V[2][i] = ZnN[i];
    }



    unsigned long **Winv;
    Winv = Allocmn3(N, N);
    unsigned long **S;
    S = Allocmn3(N, N);

    for (i = 0; i < N * SizeRowV; ++i) {
	Winv[0][i] = ZN[i];
	Winv[1][i] = ZN[i];
	Winv[2][i] = ZN[i];
	S[0][i] = IN[i];
	S[1][i] = IN[i];
	S[2][i] = IN[i];
    }


    unsigned long *SizeS;
    SizeS = malloc(3 * sizeof(unsigned long));
    SizeS[0] = N;
    SizeS[1] = N;
    SizeS[2] = N;

    unsigned long **Set;
    Set = malloc(3 * sizeof(unsigned long *));
    Set[0] = malloc((N + 1) * sizeof(unsigned long));
    Set[1] = malloc((N + 1) * sizeof(unsigned long));
    Set[2] = malloc((N + 1) * sizeof(unsigned long));

    for (i = 0; i < N; ++i) {
	Set[0][i] = i;
	Set[1][i] = i;
	Set[2][i] = i;
    }

    unsigned long *AVi;
    AVi = Allocmn(n, N);
    unsigned long *ViAVi;
    ViAVi = Allocmn(N, N);
    unsigned long *ViA2Vi;
    ViA2Vi = Allocmn(N, N);

    InnerProducts(m, n, N, a, V[0], AVi, ViAVi, ViA2Vi);

    unsigned long *UseSet;
    UseSet = malloc((N + 1) * sizeof(unsigned long));



    UseSet[0] = N;
    for (i = 1; i <= N; ++i) {
	UseSet[i] = i - 1;
    }



    unsigned long *Set0;
    Set0 = malloc(sizeof(unsigned long));

    unsigned long *nW_inv;
    nW_inv = Allocmn(N, N);




    unsigned long *ntaken;

    ntaken = malloc((N + 1) * sizeof(unsigned long));



    unsigned long *nSi;
    nSi = Allocmn(N, N);

    Construction_S_WinvBit(N, ViAVi, UseSet, nSi, Set0, nW_inv, ntaken);


    Refresh_ArrayBit(N, N, Winv, nW_inv);
    Refresh_ArrayBit(N, N, S, nSi);


    SizeS[2] = SizeS[1];
    SizeS[1] = SizeS[0];
    SizeS[0] = ntaken[0];

    for (i = 0; i < N + 1; ++i) {
	Set[2][i] = Set[1][i];
	Set[1][i] = Set[0][i];
	Set[0][i] = ntaken[i];
    }


    printf("Step= %d      SizeWi=  %lu\n", 0, SizeS[0]);

    unsigned long *AVi1;
    AVi1 = Allocmn(n, N);
    unsigned long *Vi1AVi1;
    Vi1AVi1 = Allocmn(N, N);
    unsigned long *Vi1A2Vi1;
    Vi1A2Vi1 = Allocmn(N, N);

    InnerProducts(m, n, N, a, V[1], AVi1, Vi1AVi1, Vi1A2Vi1);

//Desalloc 

    free(AVi1);


    unsigned long Step = 1;
    unsigned long SumRank = SizeS[0];


//Alloc



    unsigned long *V1i;
    V1i = Allocmn(n, N);
    unsigned long *D1i;
    D1i = Allocmn(N, N);
    unsigned long *ResNN;
    ResNN = Allocmn(N, N);
    unsigned long *ResNN2;
    ResNN2 = Allocmn(N, N);
    unsigned long *E1i;
    E1i = Allocmn(N, N);
    unsigned long *F1i;
    F1i = Allocmn(N, N);
    unsigned long *ResNN3;
    ResNN3 = Allocmn(N, N);
    unsigned long *ResnN;
    ResnN = Allocmn(n, N);
    unsigned long *ResNn;
    ResNn = Allocmn(N, n);


//the iteration


    while (TestZero(N, N, ViAVi) == 0) {


	for (i = 1; i <= Set0[0]; ++i)
	    UseSet[i] = Set[0][i];
	UseSet[0] = Set0[0];

// for V1i
	for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
	    V1i[i] = AVi[i];
// end for V1i


//        D1i:=IN-Winv`i*(ViA2Vi*S`i*Transpose(S`i)+ViAVi);

	unsigned long *ResVN;
	ResVN = Allocmn(Set0[0], N);
	unsigned long *ResNV;
	ResNV = Allocmn(N, Set0[0]);

	SelectLinesBit(N, S[0], SizeS[0], ResVN);
	TransposeBit(SizeS[0], N, ResVN, ResNV);
	DMultBit(N, SizeS[0], N, ResNV, ResVN, ResNN);


// for V1i


	DMultBit(n, N, N, V1i, ResNN, ResnN);
	for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
	    V1i[i] = ResnN[i];

// end for V1i



	DMultBit(N, N, N, ViA2Vi, ResNN, ResNN2);
	DMatrixSumBit(N, N, ResNN2, ViAVi, ResNN);
	DMultBit(N, N, N, Winv[0], ResNN, ResNN2);
	DMatrixSumBit(N, N, IN, ResNN2, D1i);



//      E1i:=-Winv`i1*ViAVi*S`i*Transpose(S`i);

	SelectLinesBit(N, S[0], SizeS[0], ResVN);
	TransposeBit(SizeS[0], N, ResVN, ResNV);
	DMultBit(N, SizeS[0], N, ResNV, ResVN, ResNN);
	DMultBit(N, N, N, ViAVi, ResNN, ResNN2);
	DMultBit(N, N, N, Winv[1], ResNN2, E1i);



// F1i:=-Winv`i2*(IN-Vi1AVi1*Winv`i1)*(Vi1A2Vi1*S`i1*Transpose(S`i1)+Vi1AVi1)*S`i*Transpose(S`i);


	unsigned long *ResVN2;
	ResVN2 = Allocmn(SizeS[1], N);
	unsigned long *ResNV2;
	ResNV2 = Allocmn(N, SizeS[1]);



	SelectLinesBit(N, S[1], SizeS[1], ResVN2);
	TransposeBit(SizeS[1], N, ResVN2, ResNV2);
	DMultBit(N, SizeS[1], N, ResNV2, ResVN2, ResNN);
	DMultBit(N, N, N, Vi1A2Vi1, ResNN, ResNN2);
	DMatrixSumBit(N, N, ResNN2, Vi1AVi1, ResNN);
	DMultBit(N, SizeS[0], N, ResNV, ResVN, ResNN2);
	DMultBit(N, N, N, ResNN, ResNN2, ResNN3);
	DMultBit(N, N, N, Vi1AVi1, Winv[1], ResNN);
	DMatrixSumBit(N, N, IN, ResNN, ResNN2);
	DMultBit(N, N, N, ResNN2, ResNN3, ResNN);
	DMultBit(N, N, N, Winv[2], ResNN, F1i);


// V1i:=AVi*S`i*Transpose(S`i)+V`i*D1i+V`i1*E1i+V`i2*F1i;


	DMultBit(n, N, N, V[0], D1i, ResnN);
	DMatrixSumBit(n, N, V1i, ResnN, V1i);
	DMultBit(n, N, N, V[1], E1i, ResnN);
	DMatrixSumBit(n, N, V1i, ResnN, V1i);
	DMultBit(n, N, N, V[2], F1i, ResnN);
	DMatrixSumBit(n, N, V1i, ResnN, V1i);



	Refresh_ArrayBit(n, N, V, V1i);




//Desalloc   


	free(ResVN);
	free(ResVN2);
	free(ResNV2);
	free(ResNV);



	for (i = 0; i < iceildiv(N, WBITS) * N; ++i) {
	    Vi1AVi1[i] = ViAVi[i];
	    Vi1A2Vi1[i] = ViA2Vi[i];
	}


	InnerProducts(m, n, N, a, V1i, AVi, ViAVi, ViA2Vi);


	Construction_S_WinvBit(N, ViAVi, UseSet, nSi, Set0, nW_inv, ntaken);

	Refresh_ArrayBit(N, N, Winv, nW_inv);
	Refresh_ArrayBit(N, N, S, nSi);

	SizeS[2] = SizeS[1];
	SizeS[1] = SizeS[0];
	SizeS[0] = ntaken[0];

	SumRank += SizeS[0];

	Refresh_ArrayBit(N + 1, 1, Set, ntaken);

	printf("Step= %lu      SizeWi=  %lu    SumRank= %lu\n", Step,
	       SizeS[0], SumRank);

	Step += 1;

	if (SumRank > n) {
	    printf("Failed\n");
	    break;
	}

    }

    for (i = 0; i < iceildiv(N, WBITS) * n; ++i)
	X[i] = V[0][i];

//Desalloc 


    free(ResNn);
    free(ResNN);
    free(ResNN2);
    free(ResnN);
    free(ResNN3);
    free(V1i);
    free(AVi);
    free(ViAVi);
    free(ViA2Vi);
    free(D1i);
    free(E1i);
    free(F1i);

    free(IN);
    free(ZN);
    free(ZnN);


    free(nW_inv);
    free(ntaken);
    free(nSi);

    free(SizeS);
    free(UseSet);
    free(Set0);

    free(Vi1AVi1);

    free(Vi1A2Vi1);

    FreeAllocmn3(V);
    FreeAllocmn3(Winv);
    FreeAllocmn3(Set);
    FreeAllocmn3(S);

    return 0;

}






////////////////////////////////////////
//
//
//   Some Preliminary functions to older versions
//
////////////////////////////////////////


void displayMatrixToto(const unsigned long *matrix,
		       unsigned long m, unsigned long n, char c)
{
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);

    printf("%c:=Matrix(GF(2),%ld,%ld,[", c, m, n);

//printf("Matrix(GF(2),%ld,%ld,[",m,n);
    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    if (i != 0 || j != 0)
		printf(", ");
	    printf("%d", bit);
	}
    }
    printf("]);\n");

}

void displayMatrixOld(unsigned long **M,
		      unsigned long m, unsigned long n, char c)
{
    unsigned long i, k;
    printf("%c:=Matrix(GF(2),%lu,%lu,[", c, m, n);
    for (k = 0; k < m; ++k) {

	for (i = 0; i < n; ++i) {
	    if ((i != n - 1) || (k != m - 1))
		printf("%lu, ", M[k][i]);
	    else
		printf("%lu ", M[k][i]);
	}
    }
    printf("]);\n");

}


unsigned long Refresh_Array(unsigned long m, unsigned long n,
			    unsigned long ***A, unsigned long **a)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i)
	for (j = 0; j < n; ++j) {
	    A[2][i][j] = A[1][i][j];
	    A[1][i][j] = A[0][i][j];
	    A[0][i][j] = a[i][j];
	}
    return 0;
}


unsigned long Mult_TV_A_V(unsigned long n, unsigned long N, unsigned long **v,
			  unsigned long **a, unsigned long **Res2)
{
    unsigned long i, j;

    unsigned long **ResnN;

    ResnN = Alloc2(n, N);

    unsigned long **ResNn;

    ResNn = Alloc2(N, n);

    SMultDmatrix(N, a, v, ResnN);
    Transpose(n, N, v, ResNn);
    DMult(N, n, N, ResNn, ResnN, Res2);


// Desalloc

    FreeAlloc2(ResnN, n);
    FreeAlloc2(ResNn, N);


    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j)
	    if (Res2[i][j] != 0)
		return 0;
    return 1;
}



unsigned long Mult_TV_A_VBit(unsigned long n, unsigned long N,
			     unsigned long *v, unsigned long *a,
			     unsigned long *Res2)
{
    unsigned long i;

    unsigned long *ResnN;

    ResnN = Allocmn(n, N);

    unsigned long *ResNn;

    ResNn = Allocmn(N, n);

    SMultDmatrixBitNew(n, N, a, v, ResnN);
    TransposeBit(n, N, v, ResNn);
    DMultBit(N, n, N, ResNn, ResnN, Res2);


    free(ResnN);
    free(ResNn);

    for (i = 0; i < N; ++i)
	if (Res2[i] != 0)
	    return 0;
    return 1;
}








/////////////////////////////////////////////////////////
//
// From now on older versions of Lanczos
//
/////////////////////////////////////////////////////////



unsigned long bigdeal(unsigned long n, unsigned long N, unsigned long **a,
		      unsigned long **v, unsigned long **X)
{
    unsigned long i, j;

    //initial conditions data prodution


    unsigned long **IN;

    IN = Alloc2(N, N);

    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j) {
	    if (i == j) {
		IN[i][j] = 1;
	    } else {
		IN[i][j] = 0;
	    }
	}


    unsigned long **ZN;

    ZN = Alloc2(N, N);

    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j) {
	    ZN[i][j] = 0;
	}


    unsigned long **ZnN;

    ZnN = Alloc2(n, N);

    for (i = 0; i < n; ++i)
	for (j = 0; j < N; ++j) {
	    ZnN[i][j] = 0;
	}

    unsigned long ***V;

    V = Alloc3(n, N);

    unsigned long ***W;

    W = Alloc3(n, N);

    for (i = 0; i < n; ++i)
	for (j = 0; j < N; ++j) {
	    V[0][i][j] = v[i][j];
	    V[1][i][j] = ZnN[i][j], V[2][i][j] = ZnN[i][j];
	    W[0][i][j] = ZnN[i][j], W[1][i][j] =
		ZnN[i][j], W[2][i][j] = ZnN[i][j];
	}
    unsigned long ***Winv;

    Winv = Alloc3(N, N);

    unsigned long ***S;

    S = Alloc3(N, N);

    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j) {
	    Winv[0][i][j] = ZN[i][j];
	    Winv[1][i][j] = ZN[i][j], Winv[2][i][j] = ZN[i][j];
	    S[0][i][j] = IN[i][j], S[1][i][j] =
		IN[i][j], S[2][i][j] = IN[i][j];
	}


    unsigned long SizeS[3][1];

    SizeS[0][0] = N;
    SizeS[1][0] = N;
    SizeS[2][0] = N;



    unsigned long ***Set;

    Set = Alloc3(N + 1, 1);

    for (i = 0; i < N; ++i) {
	Set[0][i][0] = i;
	Set[1][i][0] = i;
	Set[2][i][0] = i;
    }


    unsigned long **ResnN;

    ResnN = Alloc2(n, N);


    SMultDmatrix(N, a, v, ResnN);


    unsigned long **ResNn;

    ResNn = Alloc2(N, n);


    Transpose(n, N, v, ResNn);

    unsigned long **T;

    T = Alloc2(N, N);


    DMult(N, n, N, ResNn, ResnN, T);


//Desalloc

    FreeAlloc2(ResNn, N);
    FreeAlloc2(ResnN, n);




    unsigned long UseSet[N + 1];

    for (i = 1; i <= N; ++i)
	UseSet[i] = i - 1;

    UseSet[0] = N;


    unsigned long Set0[1];

    unsigned long **nW_inv;

    nW_inv = Alloc2(N, N);

    unsigned long **ntaken;

    ntaken = Alloc2(N + 1, 1);

    unsigned long **nSi;

    nSi = Alloc2(N, N);


    Construction_S_Winv(N, T, UseSet, nSi, Set0, nW_inv, ntaken);



    Refresh_Array(N, N, Winv, nW_inv);

    Refresh_Array(N, N, S, nSi);



    SizeS[2][0] = SizeS[1][0];
    SizeS[1][0] = SizeS[0][0];
    SizeS[0][0] = ntaken[0][0];


    Refresh_Array(N + 1, 1, Set, ntaken);



    unsigned long **Wi;

    Wi = Alloc2(n, SizeS[0][0]);


    unsigned long **ResVN;

    ResVN = Alloc2(SizeS[0][0], N);

    unsigned long **ResNV;

    ResNV = Alloc2(N, SizeS[0][0]);

    SelectLines(N, S[0], SizeS[0][0], ResVN);
    Transpose(SizeS[0][0], N, ResVN, ResNV);


    DMult(n, N, SizeS[0][0], V[0], ResNV, Wi);

    Refresh_Array(n, SizeS[0][0], W, Wi);



//Desalloc

    FreeAlloc2(Wi, n);
    FreeAlloc2(ResNV, N);
    FreeAlloc2(ResVN, SizeS[0][0]);






    //printf("Step= %u      SizeWi=  %u\n", 0, SizeS[0][0]);


    unsigned long Step = 1;
    unsigned long SumRank = SizeS[0][0];

    unsigned long **Res2;

    Res2 = Alloc2(N, N);


    unsigned long **try;

    try = Alloc2(n, N);


//the iteration


    while (Mult_TV_A_V(n, N, V[0], a, Res2) == 0) {

	unsigned long UseSet[N + 1];
	for (i = 1; i <= Set0[0]; ++i)
	    UseSet[i] = Set[0][i][0];
	UseSet[0] = Set0[0];



	unsigned long **Res4;

	Res4 = Alloc2(n, N);

	SMultDmatrix(N, a, V[0], Res4);

	unsigned long **ResVN;

	ResVN = Alloc2(SizeS[0][0], N);

	unsigned long **ResNV;

	ResNV = Alloc2(N, SizeS[0][0]);

	unsigned long **ResnV;

	ResnV = Alloc2(n, SizeS[0][0]);

	unsigned long **ResnN;

	ResnN = Alloc2(n, N);

	SelectLines(N, S[0], SizeS[0][0], ResVN);
	Transpose(SizeS[0][0], N, ResVN, ResNV);
	DMult(n, N, SizeS[0][0], Res4, ResNV, ResnV);
	DMult(n, SizeS[0][0], N, ResnV, ResVN, ResnN);
	DSum(n, N, ResnN, V[0], try);


//Desalloc

	FreeAlloc2(ResnV, n);
	FreeAlloc2(ResnN, n);
	FreeAlloc2(ResNV, N);
	FreeAlloc2(ResVN, SizeS[0][0]);



	unsigned long **V1i;

	V1i = Alloc2(n, N);


	for (i = 0; i < n; ++i)
	    for (j = 0; j < N; ++j)
		V1i[i][j] = try[i][j];


	SMultDmatrix(N, a, try, Res4);

	unsigned long **ResNN;

	ResNN = Alloc2(N, N);

	ResnN = Alloc2(n, N);

	unsigned long **Res2nN;

	Res2nN = Alloc2(n, N);

	unsigned long **ResNn;

	ResNn = Alloc2(N, n);

	for (i = 0; i < 3; ++i) {
	    DMult(n, N, N, V[i], Winv[i], ResnN);

	    Transpose(n, N, V[i], ResNn);

	    DMult(N, n, N, ResNn, Res4, ResNN);

	    DMult(n, N, N, ResnN, ResNN, Res2nN);

	    DSum(n, N, V1i, Res2nN, V1i);

	}

//Desalloc

	FreeAlloc2(ResnN, n);
	FreeAlloc2(Res2nN, n);
	FreeAlloc2(ResNN, N);
	FreeAlloc2(ResNn, N);




	Refresh_Array(n, N, V, V1i);

	unsigned long **Res3;

	Res3 = Alloc2(N, N);

	Mult_TV_A_V(n, N, V1i, a, Res3);

	Construction_S_Winv(N, Res3, UseSet, nSi, Set0, nW_inv, ntaken);

	FreeAlloc2(Res3, N);

	Refresh_Array(N, N, Winv, nW_inv);
	Refresh_Array(N, N, S, nSi);
	SizeS[2][0] = SizeS[1][0];
	SizeS[1][0] = SizeS[0][0];
	SizeS[0][0] = ntaken[0][0];

	SumRank += SizeS[0][0];

	Refresh_Array(N + 1, 1, Set, ntaken);




	unsigned long **Wi;

	Wi = Alloc2(n, SizeS[0][0]);

	ResVN = Alloc2(SizeS[0][0], N);

	ResNV = Alloc2(N, SizeS[0][0]);

	SelectLines(N, S[0], SizeS[0][0], ResVN);
	Transpose(SizeS[0][0], N, ResVN, ResNV);

	DMult(n, N, SizeS[0][0], V[0], ResNV, Wi);

	Refresh_Array(n, SizeS[0][0], W, Wi);

//Desalloc

	FreeAlloc2(Wi, n);
	FreeAlloc2(ResNV, N);
	FreeAlloc2(ResVN, SizeS[0][0]);



	printf("Step= %lu      SizeWi=  %lu    SumRank= %lu\n", Step,
	       SizeS[0][0], SumRank);

	Step += 1;


//Desaloc

	FreeAlloc2(V1i, n);
	FreeAlloc2(Res4, n);


    }

    for (i = 0; i < n; ++i)
	for (j = 0; j < N; ++j)
	    X[i][j] = V[0][i][j];

//Desalloc

    FreeAlloc2(IN, N);
    FreeAlloc2(ZN, N);
    FreeAlloc2(ZnN, n);

    FreeAlloc3(V, n);
    FreeAlloc3(W, n);
    FreeAlloc3(Winv, N);
    FreeAlloc3(Set, N + 1);
    FreeAlloc3(S, N);

    FreeAlloc2(Res2, N);
    FreeAlloc2(T, N);
    FreeAlloc2(nW_inv, N);
    FreeAlloc2(ntaken, N + 1);
    FreeAlloc2(nSi, N);
    FreeAlloc2(try, n);


    return 0;

}



























////////////////////////////////////////////////////
//
//
//  The Lanczos Algorithm  (Faster Version with a level two recursion
//
//
/////////////////////////////////////////////////////



unsigned long bigdeal_Fast(unsigned long n, unsigned long N,
			   unsigned long **a, unsigned long **v,
			   unsigned long **X)
{
    unsigned long i, j;

    //initial conditions data prodution


    unsigned long **IN;
    IN = Alloc2(N, N);

    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j) {
	    if (i == j) {
		IN[i][j] = 1;
	    } else {
		IN[i][j] = 0;
	    }
	}


    unsigned long **ZN;
    ZN = Alloc2(N, N);

    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j) {
	    ZN[i][j] = 0;
	}


    unsigned long **ZnN;
    ZnN = Alloc2(n, N);

    for (i = 0; i < n; ++i)
	for (j = 0; j < N; ++j) {
	    ZnN[i][j] = 0;
	}

    unsigned long ***V;
    V = Alloc3(n, N);
    unsigned long ***W;
    W = Alloc3(n, N);

    for (i = 0; i < n; ++i)
	for (j = 0; j < N; ++j) {
	    V[0][i][j] = v[i][j];
	    V[1][i][j] = ZnN[i][j], V[2][i][j] = ZnN[i][j];
	    W[0][i][j] = ZnN[i][j], W[1][i][j] =
		ZnN[i][j], W[2][i][j] = ZnN[i][j];
	}
    unsigned long ***Winv;
    Winv = Alloc3(N, N);
    unsigned long ***S;
    S = Alloc3(N, N);

    for (i = 0; i < N; ++i)
	for (j = 0; j < N; ++j) {
	    Winv[0][i][j] = ZN[i][j];
	    Winv[1][i][j] = ZN[i][j], Winv[2][i][j] = ZN[i][j];
	    S[0][i][j] = IN[i][j], S[1][i][j] =
		IN[i][j], S[2][i][j] = IN[i][j];
	}


    unsigned long SizeS[3][1];

    SizeS[0][0] = N;
    SizeS[1][0] = N;
    SizeS[2][0] = N;



    unsigned long ***Set;
    Set = Alloc3(N + 1, 1);

    for (i = 0; i < N; ++i) {
	Set[0][i][0] = i;
	Set[1][i][0] = i;
	Set[2][i][0] = i;
    }


    unsigned long **ResnN;
    ResnN = Alloc2(n, N);
    unsigned long **ResNn;
    ResNn = Alloc2(N, n);
    unsigned long **T;
    T = Alloc2(N, N);

    SMultDmatrix(N, a, v, ResnN);
    Transpose(n, N, v, ResNn);
    DMult(N, n, N, ResNn, ResnN, T);


//Desalloc

    // FreeAlloc2(ResNn,N);
    //  FreeAlloc2(ResnN,n);




    unsigned long UseSet[N + 1];

    for (i = 1; i <= N; ++i)
	UseSet[i] = i - 1;

    UseSet[0] = N;


    unsigned long Set0[1];
    unsigned long **nW_inv;
    nW_inv = Alloc2(N, N);
    unsigned long **ntaken;
    ntaken = Alloc2(N + 1, 1);
    unsigned long **nSi;
    nSi = Alloc2(N, N);


    Construction_S_Winv(N, T, UseSet, nSi, Set0, nW_inv, ntaken);



    Refresh_Array(N, N, Winv, nW_inv);
    Refresh_Array(N, N, S, nSi);



    SizeS[2][0] = SizeS[1][0];
    SizeS[1][0] = SizeS[0][0];
    SizeS[0][0] = ntaken[0][0];


    Refresh_Array(N + 1, 1, Set, ntaken);


    unsigned long **Wi;
    Wi = Alloc2(n, SizeS[0][0]);
    unsigned long **ResVN;
    ResVN = Alloc2(SizeS[0][0], N);
    unsigned long **ResNV;
    ResNV = Alloc2(N, SizeS[0][0]);

    SelectLines(N, S[0], SizeS[0][0], ResVN);
    Transpose(SizeS[0][0], N, ResVN, ResNV);


    DMult(n, N, SizeS[0][0], V[0], ResNV, Wi);

    Refresh_Array(n, SizeS[0][0], W, Wi);



//Desalloc

    FreeAlloc2(Wi, n);
    FreeAlloc2(ResNV, N);
    FreeAlloc2(ResVN, SizeS[0][0]);



    printf("Step= %d      SizeWi=  %lu\n", 0, SizeS[0][0]);


    unsigned long **AVi1;
    AVi1 = Alloc2(n, N);
    unsigned long **Vi1AVi1;
    Vi1AVi1 = Alloc2(N, N);
    unsigned long **Vi1A2Vi1;
    Vi1A2Vi1 = Alloc2(N, N);


    SMultDmatrix(N, a, V[1], AVi1);
    Transpose(n, N, V[1], ResNn);
    DMult(N, n, N, ResNn, AVi1, Vi1AVi1);
    Transpose(n, N, AVi1, ResNn);
    DMult(N, n, N, ResNn, AVi1, Vi1A2Vi1);

//Desalloc

    FreeAlloc2(AVi1, n);




    unsigned long Step = 1;
    unsigned long SumRank = SizeS[0][0];



//Alloc

    unsigned long **Res2;
    Res2 = Alloc2(N, N);
    unsigned long **AVi;
    AVi = Alloc2(n, N);
    unsigned long **ViAVi;
    ViAVi = Alloc2(N, N);
    unsigned long **ViA2Vi;
    ViA2Vi = Alloc2(N, N);
    unsigned long **V1i;
    V1i = Alloc2(n, N);
    unsigned long **D1i;
    D1i = Alloc2(N, N);
    unsigned long **ResNN;
    ResNN = Alloc2(N, N);
    unsigned long **ResNN2;
    ResNN2 = Alloc2(N, N);
    unsigned long **E1i;
    E1i = Alloc2(N, N);
    unsigned long **F1i;
    F1i = Alloc2(N, N);
    unsigned long **ResNN3;
    ResNN3 = Alloc2(N, N);
    unsigned long **Res3;
    Res3 = Alloc2(N, N);







//the iteration


    while (Mult_TV_A_V(n, N, V[0], a, Res2) == 0) {

	unsigned long UseSet[N + 1];
	for (i = 1; i <= Set0[0]; ++i)
	    UseSet[i] = Set[0][i][0];
	UseSet[0] = Set0[0];



	SMultDmatrix(N, a, V[0], AVi);

	Transpose(n, N, V[0], ResNn);
	DMult(N, n, N, ResNn, AVi, ViAVi);

	Transpose(n, N, AVi, ResNn);
	DMult(N, n, N, ResNn, AVi, ViA2Vi);

	for (i = 0; i < n; ++i)
	    for (j = 0; j < N; ++j) {
		V1i[i][j] = AVi[i][j];
	    }

//   Verificado






//        D1i:=IN-Winv`i*(ViA2Vi*S`i*Transpose(S`i)+ViAVi);

	unsigned long **ResVN;
	ResVN = Alloc2(Set0[0], N);
	unsigned long **ResNV;
	ResNV = Alloc2(N, Set0[0]);

	SelectLines(N, S[0], SizeS[0][0], ResVN);

	Transpose(SizeS[0][0], N, ResVN, ResNV);
	DMult(N, SizeS[0][0], N, ResNV, ResVN, ResNN);




// for V1i

	DMult(n, N, N, V1i, ResNN, ResnN);



	for (i = 0; i < n; ++i)
	    for (j = 0; j < N; ++j)
		V1i[i][j] = ResnN[i][j];


// end for V1i


	DMult(N, N, N, ViA2Vi, ResNN, ResNN2);



	DSum(N, N, ResNN2, ViAVi, ResNN);
	DMult(N, N, N, Winv[0], ResNN, ResNN2);
	DSum(N, N, IN, ResNN2, D1i);






//      E1i:=-Winv`i1*ViAVi*S`i*Transpose(S`i);

	SelectLines(N, S[0], SizeS[0][0], ResVN);
	Transpose(SizeS[0][0], N, ResVN, ResNV);
	DMult(N, SizeS[0][0], N, ResNV, ResVN, ResNN);
	DMult(N, N, N, ViAVi, ResNN, ResNN2);
	DMult(N, N, N, Winv[1], ResNN2, E1i);




// F1i:=-Winv`i2*(IN-Vi1AVi1*Winv`i1)*(Vi1A2Vi1*S`i1*Transpose(S`i1)+Vi1AVi1)*S`i*Transpose(S`i);


	unsigned long **ResVN2;
	ResVN2 = Alloc2(SizeS[1][0], N);
	unsigned long **ResNV2;
	ResNV2 = Alloc2(N, SizeS[1][0]);



	SelectLines(N, S[1], SizeS[1][0], ResVN2);
	Transpose(SizeS[1][0], N, ResVN2, ResNV2);
	DMult(N, SizeS[1][0], N, ResNV2, ResVN2, ResNN);
	DMult(N, N, N, Vi1A2Vi1, ResNN, ResNN2);
	DSum(N, N, ResNN2, Vi1AVi1, ResNN);
	DMult(N, SizeS[0][0], N, ResNV, ResVN, ResNN2);
	DMult(N, N, N, ResNN, ResNN2, ResNN3);
	DMult(N, N, N, Vi1AVi1, Winv[1], ResNN);
	DSum(N, N, IN, ResNN, ResNN2);
	DMult(N, N, N, ResNN2, ResNN3, ResNN);
	DMult(N, N, N, Winv[2], ResNN, F1i);




//displayMatrixOld(F1i,N,N,'f');
	// V1i:=AVi*S`i*Transpose(S`i)+V`i*D1i+V`i1*E1i+V`i2*F1i;


	DMult(n, N, N, V[0], D1i, ResnN);
	DSum(n, N, V1i, ResnN, V1i);
	DMult(n, N, N, V[1], E1i, ResnN);
	DSum(n, N, V1i, ResnN, V1i);
	DMult(n, N, N, V[2], F1i, ResnN);
	DSum(n, N, V1i, ResnN, V1i);


//displayMatrixOld(V1i,n,N,'f');

//      displayMatrixOld(V1i,n,N,'e');

	Refresh_Array(n, N, V, V1i);

//      displayMatrixOld(V[0],n,N,'f');

//Desalloc


	FreeAlloc2(ResVN, Set0[0]);
	FreeAlloc2(ResVN2, SizeS[1][0]);
	FreeAlloc2(ResNV2, N);
	FreeAlloc2(ResNV, N);




	for (i = 0; i < N; ++i)
	    for (j = 0; j < N; ++j) {
		Vi1AVi1[i][j] = ViAVi[i][j];
		Vi1A2Vi1[i][j] = ViA2Vi[i][j];
	    }



	Mult_TV_A_V(n, N, V1i, a, Res3);



	Construction_S_Winv(N, Res3, UseSet, nSi, Set0, nW_inv, ntaken);



	Refresh_Array(N, N, Winv, nW_inv);
	Refresh_Array(N, N, S, nSi);
	SizeS[2][0] = SizeS[1][0];
	SizeS[1][0] = SizeS[0][0];
	SizeS[0][0] = ntaken[0][0];

	SumRank += SizeS[0][0];

	Refresh_Array(N + 1, 1, Set, ntaken);

	unsigned long **Wi;
	Wi = Alloc2(n, SizeS[0][0]);
	ResVN = Alloc2(SizeS[0][0], N);
	ResNV = Alloc2(N, SizeS[0][0]);

	SelectLines(N, S[0], SizeS[0][0], ResVN);
	Transpose(SizeS[0][0], N, ResVN, ResNV);

	DMult(n, N, SizeS[0][0], V[0], ResNV, Wi);



	Refresh_Array(n, SizeS[0][0], W, Wi);

//Desalloc

	FreeAlloc2(Wi, n);
	FreeAlloc2(ResNV, N);
	FreeAlloc2(ResVN, SizeS[0][0]);



	printf("Step= %lu      SizeWi=  %lu    SumRank= %lu\n", Step,
	       SizeS[0][0], SumRank);

	Step += 1;



    }


    for (i = 0; i < n; ++i)
	for (j = 0; j < N; ++j)

	    X[i][j] = V[0][i][j];

//Desalloc

    FreeAlloc2(Res3, N);
    FreeAlloc2(ResNn, N);
    FreeAlloc2(ResNN, N);
    FreeAlloc2(ResNN2, N);
    FreeAlloc2(ResnN, n);
    FreeAlloc2(ResNN3, N);
    FreeAlloc2(V1i, n);
    FreeAlloc2(AVi, n);
    FreeAlloc2(ViAVi, N);
    FreeAlloc2(ViA2Vi, N);
    FreeAlloc2(D1i, N);
    FreeAlloc2(E1i, N);
    FreeAlloc2(F1i, N);

    FreeAlloc2(IN, N);
    FreeAlloc2(ZN, N);
    FreeAlloc2(ZnN, n);
    FreeAlloc2(Res2, N);
    FreeAlloc2(T, N);
    FreeAlloc2(nW_inv, N);
    FreeAlloc2(ntaken, N + 1);
    FreeAlloc2(nSi, N);

    FreeAlloc2(Vi1AVi1, N);
    FreeAlloc2(Vi1A2Vi1, N);

    FreeAlloc3(V, n);
    FreeAlloc3(W, n);
    FreeAlloc3(Winv, N);
    FreeAlloc3(Set, N + 1);
    FreeAlloc3(S, N);

    return 0;

}
