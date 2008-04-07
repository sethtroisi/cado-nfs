#include <iomanip>
#include "bitstring.hpp"
#include "manu.h"
#include "utils.h"

using namespace std;
unsigned int read_hexstring(std::istream& is, unsigned long * ptr, unsigned int n)
{
    unsigned int i = 0;
    unsigned int ir = 0;
    unsigned long x = 0UL;
    unsigned long v;
    is >> ws;
    for( ; i < n ; ) {
        char c;
        c = is.peek();
        if (c >= '0' && c <= '9') v = c - '0';
        else if (c >= 'a' && c <= 'f') v = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v = c - 'A' + 10;
        else {
            return i;
        }
        is.ignore(1);

        x |= v << ir;
        i += 4;
        ir += 4;
        if (ir == ULONG_BITS) {
            ir = 0;
            *ptr++ = x;
            x = 0UL;
        }
    }
    if (i == n)
        return n;

    /* now n < i < n + 4. According to the specification, we check that
     * the top i-n bits are zero, and return n.
     */
    v >>= n+4-i;
    BUG_ON(v);
    return n;
}

unsigned int write_hexstring(std::ostream& o, const unsigned long * ptr, unsigned int n)
{
    const char hxits[] = "0123456789abcdef";
    unsigned int i = 0;
    for( ; i < n ; ) {
        unsigned long v = *ptr++;
        for(unsigned int ir = 0 ; ir < ULONG_BITS && i < n ; ) {
            o << hxits[v&0x0f];
            if (!o) return i;
            v >>= 4;
            ir += 4;
            i += 4;
        }
    }
    return n;
}
