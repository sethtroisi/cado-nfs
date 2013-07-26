/* This source file is our test case for I/O odd adress for an int. If this can be
 * compiled, it means that ODD_ADDRESS_IO_INT will be usable in sieve/las.
 */

int main() {
  unsigned long long i = 0x0123456789ABCDEF;
  unsigned int j = *(unsigned int *) (((unsigned char *) &i) + 1);
  /* Little endian & big endian are possible */
  return (j == 0x6789ABCD || j == 0x23456789) ? 0 : 1;
}

