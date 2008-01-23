#include <stdio.h>
#include <stdlib.h>

#define PRINTVAL(x)	do {		\
	if (!first) printf(", ");	\
		first=0;		\
		printf("%d",x);		\
	} while (0)
	
#define MAGMA

#ifdef MAGMA
#define START_MATRIX	started_matrix=0
#define START_LINE	do {		\
	if (!started_matrix) {		\
		printf("[ ");		\
		started_matrix=1;	\
	} else {			\
		printf(", ");		\
	}				\
	first=1;			\
	} while (0)

#define END_LINE
#define END_MATRIX	printf(" ];\n");
#else
#define START_MATRIX
#define START_LINE	do {		\
	first=1;			\
	printf("[ ");			\
	} while (0)
#define END_LINE	printf(" ]\n");
#define END_MATRIX	printf("\n");
#endif
		
int main(int argc, char *argv[])
{
	char buffer[100];
	FILE *f;
	int len;
	unsigned int idx;
	int coeff;
	int ucoeff;
	int i, ii, j, k, first, started_matrix=0;
	int bank;
	char filename[80];
	int m, n, size;

	if (argc != 5) {
		fprintf(stderr,"Usage : dumpmat <bank#> <N> <m> <n>\n");
		exit(1);
	}

	bank = atoi(argv[1]);
	size = atoi(argv[2]);
	m  = atoi(argv[3]);
	n  = atoi(argv[4]);

	sprintf(filename,"run/input/matrix.%d",bank);
	
	printf("// MATRIX : \n// ");

	f = fopen(filename,"r");

	fgets(buffer,100,f);
	printf(buffer);
	printf("\n\nMlist:=");
	START_MATRIX;
	for(;!feof(f);) {
		fread(&len,sizeof(int),1,f);
		if (feof(f)) break;
		ii=0;
		START_LINE;
		for(i=0;i<len;i++) {
			fread(&idx, sizeof(int),1,f);
			fread(&coeff, sizeof(int),1,f);
			for(;ii<idx;ii++)
				PRINTVAL(0);
			PRINTVAL(coeff);
			ii++;
		}

		for(;ii<size;ii++)
			PRINTVAL(0);
		END_LINE;
	}
	END_MATRIX;

	fclose(f);

	sprintf(filename,"run/input/vector.%d",bank);
	printf("\n\n// VECTOR : \n// ");

	f = fopen(filename,"r");

	fgets(buffer,100,f);
	printf(buffer);
	printf("\n\nVlist:=");
	START_MATRIX;
	for(;!feof(f);) {
		fread(&coeff, sizeof(int),1,f);
		if (feof(f)) break;
		START_LINE;
		PRINTVAL(coeff);
		END_LINE;
	}
	END_MATRIX;
	fclose(f);

#ifdef MAGMA
	printf("Zlist:=[[] : i in [1..%d]];\n",n);
#endif
	for(i=0;i<n;i++) {
		sprintf(filename,"run/output/%d/Z<%02d>",bank,i);
		f = fopen(filename,"r");
		if (f==NULL) continue;
		printf("\n\n// VECTOR Z%d : \n",i);
		printf("Zlist[%d]:=",i+1);
		START_MATRIX;
		for(;!feof(f);) {
			fread(&coeff, sizeof(int),1,f);
			if (feof(f)) break;
			START_LINE;
			PRINTVAL(coeff);
			END_LINE;
		}
		END_MATRIX;
		fclose(f);
	}

#ifdef MAGMA
	printf("Ylist:=[[] : i in [1..%d]];\n",n);
#endif
	for(i=0;i<n;i++) {
		sprintf(filename,"run/output/%d/Y<%02d>",bank,i);
		f = fopen(filename,"r");
		if (f==NULL) continue;
		printf("\n\n// VECTOR Y%d : \n",i);
		printf("Ylist[%d]:=",i+1);
		START_MATRIX;
		for(;!feof(f);) {
			fread(&coeff, sizeof(int),1,f);
			if (feof(f)) break;
			START_LINE;
			PRINTVAL(coeff);
			END_LINE;
		}
		END_MATRIX;
		fclose(f);
	}

#ifdef MAGMA
	printf("
origN:=%d;
N:=%d;
m:=%d;
n:=%d;
p:=997;
did_inits:=1;

load\"%s/BlockWiedemann/Magma/block-wiedemann.mg\";
CX:=rec<CONTEXT|>;
IT:=rec<ITERATOR|>;

Blist:=(&cat [[Vlist[i + 1]] cat [Mlist[1 + i*origN+j] : j in [0..origN-1]] : i in
[0..origN-1]]) cat [0 : i in [0..N-1]];
B:=KMN!Blist;
CX`x:=KMNm!0;
for i:=1 to m do CX`x[i,i]:=1; end for;
 
CX`y:=KMNn!(Transpose(KMatrixSpace(K,n,N)!(&cat [Ylist[i] : i in [1..n]])));
CX`z:=KMNn!(Transpose(KMatrixSpace(K,n,N)!(&cat [Zlist[i] : i in [1..n]])));

for i:=1 to N do
	blah:=Integers()!B[i,1];

	if (p - blah) le blah then
		for j:=1 to N do
			B[i,j]:=-B[i,j];
		end for;
	end if;
end for;
		
",size,size+1,m,n,getenv("HOME"));
	
#endif

	printf("CX`a:=KPMnm!0;\n");
	for(i=0;i<m;i++) for(j=0;j<n;j++) {
		sprintf(filename,"run/output/%d/A<%02d,%02d>.000",bank,i,j);
		/* On se limite à 000 parce que sinon, c'est un peu dur à
		 * voir à la main de toute façon. */
		f=fopen(filename,"r");
		if (f==NULL) continue;
		/* XXX  : On prend la transposée !!! XXX */
		printf("CX`a[%d,%d]:=",j+1,i+1);
		for(k=0;;k++) {
			fread(&ucoeff, sizeof(int),1,f);
			if (feof(f)) break;
			if (k) printf(" + ");
			printf("%u * X^%d", ucoeff, k);
		}
		printf(";\n");
		fclose(f);
	}

	return 0;
}
