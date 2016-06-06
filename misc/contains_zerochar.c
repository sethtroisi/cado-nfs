#include <stdio.h>
#include <stdlib.h>

int main() {
    int c;
    int pos = 0;
    while ((c=getchar()) != EOF) {
        if (c == '\0') {
            printf("Got a zero char at position %d.\n", pos);
            return 1;
        }
        pos++;
    }
    return 0;
}


