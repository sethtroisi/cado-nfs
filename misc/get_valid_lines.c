#include <stdio.h>
#include <string.h>
#include <assert.h>

#define N 1024

int
main()
{
  int c, ok;
  char s[N], *r;
  unsigned long lines = 0, valid = 0, invalid = 0, l;

  while ((c = getchar ()) != EOF)
    {
      ungetc (c, stdin);
      r = fgets (s, N, stdin);
      if (r == NULL)
        break;
      lines ++;
      l = strlen (s);
      assert (l < N - 1);
      /* check s[l-1] is '\n' */
      assert (s[l-1] == '\n');
      s[l-1] = '\0';

      r = s;
      ok = 0;

      /* check a */
      if (r[0] == '-')
        r++;
      if (!isdigit (r[0]))
        goto skip;
      do r++; while (isdigit (r[0]));

      /* check ',' */
      if (r[0] != ',')
        goto skip;
      r++; /* read comma */

      /* check b */
      if (!isdigit (r[0]))
        goto skip;
      do r++; while (isdigit (r[0]));

      /* check ':' */
      if (r[0] != ':')
        goto skip;
      r++; /* read ':' */

      /* check at least one rational prime */
      if (!isxdigit (r[0]))
        goto skip;
      do r++; while (isxdigit (r[0]));
      
      while (r[0] == ',')
        {
          r++; /* skip ',' */
          /* read another rational prime */
          if (!isxdigit (r[0]))
            goto skip;
          do r++; while (isxdigit (r[0]));
        }

      /* check ':' */
      if (r[0] != ':')
        goto skip;
      r++;

      /* check at least one algebraic prime */
      if (!isxdigit (r[0]))
        goto skip;
      do r++; while (isxdigit (r[0]));
      
      while (r[0] == ',')
        {
          r++; /* skip ',' */
          /* read another algebraic prime */
          if (!isxdigit (r[0]))
            goto skip;
          do r++; while (isxdigit (r[0]));
        }

      /* check we have read the whole line */
      if (r != s + strlen (s))
        goto skip;

      ok = 1; /* valid line */

    skip:
      if (ok)
        {
          printf ("%s\n", s);
          valid ++;
        }
      else
        {
          fprintf (stderr, "Skip invalid line: %s\n", s);
          invalid ++;
        }
    }
  fprintf (stderr, "Read %lu lines, %lu valid, %lu invalid\n",
           lines, valid, invalid);
  return 0;
}
