#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    FILE *in, *out;
    unsigned long x;
    int i, nbits = 8 * sizeof(unsigned long), nb, b;

    out = fopen(argv[1], "w");
    for(i = 2; i < argc; i++){
	fprintf(stderr, "Operating on file %s\n", argv[i]);
	// format is one dependency per file
	in = fopen(argv[i], "r");
	// one bit per line
	x = 0;
	nb = 0;
	while(fscanf(in, "%x", &b) != EOF){
	    if(b)
		x |= (1UL) << nb;
	    nb++;
	    if(nb == nbits){
		fprintf(out, "%lx ", x);
		x = 0;
		nb = 0;
	    }
	}
	if(nb)
	    fprintf(out, "%lx", x);
	fprintf(out, "\n");
	fclose(in);
    }
    fclose(out);
    return 0;
}
