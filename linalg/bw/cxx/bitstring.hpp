#ifndef BITSTRING_HPP_
#define BITSTRING_HPP_

#include <istream>
#include <ostream>

/* These functions read bitstrings in hex format.
 *
 * The following are guaranteed.
 * - no spaces are printed.
 * - read_hexstring() returns the number of bits represented in the
 *   string, rounded up to a multiple of four -- or possibly not rounded
 *   if the length argument wasn't rounded either. In the latter case,
 *   trailing bits must be zero, or failure will occur.
 * - read_hexstring(o, p, n1*w) followed by read_hexstring(o,p+n1,n2*w)
 *   has the exact same effect as read_hexstring(o, p, (n1+n2)*w). Here w
 *   is ULONG_BITS.
 *
 * The consequence of these requirements (which are truly useful for the
 * bw code) is that the data direction might seem odd.
 *
 * ab987531 gives the unsigned long 0x135789ba
 *
 * So the first nibble is the lowest nibble of the resulting ulong
 * string. It's quite disturbing, yes. It's extreme-little-endian, in
 * some sense (the mere fact of writing 0x135789ba is a convention which
 * could be debated).
 */
unsigned int read_hexstring(std::istream&, unsigned long *, unsigned int);
unsigned int write_hexstring(std::ostream&, const unsigned long *, unsigned int);
unsigned int read_hexstring_u64(std::istream&, uint64_t *, unsigned int);
unsigned int write_hexstring_u64(std::ostream&, const uint64_t *, unsigned int);
unsigned int read_hexstring(FILE * f, unsigned long *, unsigned int);
unsigned int write_hexstring(FILE * f, const unsigned long *, unsigned int);

#endif	/* BITSTRING_HPP_ */
