/* Do not check that p is really a bad ideal, must call renumber_is_bad before */
static inline void
handle_bad_ideals (MAYBE_UNUSED int *exp_above, int64_t a, uint64_t b,
                   unsigned long p, MAYBE_UNUSED int e)
{
#if 0
  /* handle bad ideals for the following polynomial:
      c4: 120
      c3: 467713
      c2: -2284154258
      c1: 57498079368420
      c0: -21758746190334120
    Bad ideals are:
      # I1: [ [2, 0, 0, 0], [0, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 1] ]
      # I2: [ [2, 0, 0, 0], [0, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1] ]
      2,0:1: 2
      # I3: [ [3, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1] ]
      # I4: [ [3, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [2, 0, 0, 1] ]
      3,0:1: 2
    Used for a p59 NFS-DL
  */
  if (p == 2)
  {
    int amod8 = a%8;
    if (amod8 < 0)
      amod8 += 8;
    int bmod8 = b%8;
    ASSERT_ALWAYS(amod8 == 0 || amod8 == 2 || amod8 == 4 || amod8 == 6);
    ASSERT_ALWAYS(bmod8 == 1 || bmod8 == 3 || bmod8 == 5 || bmod8 == 7);
    int v = (amod8*bmod8) % 8;
    if (v == 0 || v == 4)
    {
      ASSERT_ALWAYS(e == 3);
      exp_above[0] = 1;
      exp_above[1] = 2;
    }
    else if (v == 2)
    {
      ASSERT_ALWAYS(e == 5);
      exp_above[0] = 2;
      exp_above[1] = 3;
    }
    else if (v == 6)
    {
      ASSERT_ALWAYS(e >= 6);
      exp_above[0] = e-3;
      exp_above[1] = 3;
    }
    else
      ASSERT_ALWAYS(0);
  }
  else if (p == 3)
  {
    int amod9 = a%9;
    if (amod9 < 0)
      amod9 += 9;
    int bmod9 = b%9;
    ASSERT_ALWAYS(amod9 == 0 || amod9 == 3 || amod9 == 6);
    ASSERT_ALWAYS(bmod9 == 1 || bmod9 == 2 || bmod9 == 4 || bmod9 == 5 || bmod9 == 7 || bmod9 == 8);
    if (bmod9 == 2 || bmod9 == 4)
      bmod9 += 3;
    else if (bmod9 == 5 || bmod9 == 7)
      bmod9 -= 3;
    int v = (amod9*bmod9) % 9;
    if (v == 0)
    {
      ASSERT_ALWAYS(e == 2);
      exp_above[0] = 1;
      exp_above[1] = 1;
    }
    else if (v == 3)
    {
      ASSERT_ALWAYS(e >= 2);
      exp_above[0] = 1;
      exp_above[1] = e-1;
    }
    else if (v == 6)
    {
      ASSERT_ALWAYS(e >= 2);
      exp_above[0] = e-1;
      exp_above[1] = 1;
    }
    else
      ASSERT_ALWAYS(0);
  }
  else
    ASSERT_ALWAYS(0);

#elif 1
  /* handle bad ideals for the following polynomial:
      c5: 1919367450
      c4: -372912695938455
      c3: 28914191881484128958
      c2: 2214553172801635020339838
      c1: 16533777180271226594332227762
      c0: -1218572890737513374062358400776760
    Bad ideals are:
      2,0:1: 2
      ## Columns for (2,0):   always [ e-1, 1 ]
      # More precisely, if a/b mod 4 = 2 ->  [ 1, 1 ]
      #                    a/b mod 8 = 0 ->  [ 2, 1 ]
      #                    a/b mod 16 = 4 -> [ 3, 1 ]
      #                    a/b mod 16 = 12 -> [e-1, 1]
      5,2:1: 2
      ## Columns for (5,2):
      # a/b mod 25 in {2,7,22}  ->  [ 1,   1   ]
      # a/b mod 25 == 12        ->  [ 1,   e-1 ]
      # a/b mod 25 == 17        ->  [ e-1, 1   ]
    Used for a p155 NFS-DL
  */
  if (p == 2)
  {
    exp_above[0] = e-1;
    exp_above[1] = 1;
#ifdef DEBUG
    //hack: this compute a/b % 16 (not very efficient)
    unsigned long r = findroot (a, b, 16);
    if (r == 2 || r == 6 || r == 10 || r == 14)
      ASSERT_ALWAYS(e == 2);
    else if (r == 0 || r == 8)
      ASSERT_ALWAYS(e == 3);
    else if (r == 4)
      ASSERT_ALWAYS(e == 4);
    else if (r == 12)
      ASSERT_ALWAYS(1); //which assert here ??
    else
      ASSERT_ALWAYS(0);
#endif
  }
  else if (p == 5)
  {
    //hack: this compute a/b % 16 (not very efficient)
    unsigned long r = findroot (a, b, 25);
    if (r == 2 || r == 7 || r == 22)
    {
      exp_above[0] = 1;
      exp_above[1] = 1;
    }
    else if (r == 12)
    {
      exp_above[0] = 1;
      exp_above[1] = e-1;
    }
    else if (r == 17)
    {
      exp_above[0] = e-1;
      exp_above[1] = 1;
    }
    else
      ASSERT_ALWAYS(0);
  }
  else
    ASSERT_ALWAYS(0);
#elif 0
  /* handle bad ideals for the following polynomial:
        pol0=19B3,326,1AB,6B,0,7,1
    Bad ideals are:
      7,2:0: 2
    Used for FFS for 2^137 and 2^809
  */

  ASSERT_ALWAYS(p == 7);

  fppol64_t pol_a, pol_b;
  fppol64_t seven; // 7 is the nickname of t^2+t+1

  fppol64_set_ui_sparse (pol_a, (uint64_t) a);
  fppol64_set_ui_sparse (pol_b, b);
  fppol64_set_ui_sparse (seven, 7);

  fppol64_t Dpol_a, Dpol_b;  // derivatives of a and b
  uint64_t mask = 6148914691236517205U; // 1+t^2+t^4+...
  fppol64_div_ti(Dpol_a, pol_a, 1);
  Dpol_a[0] &= mask;
  fppol64_div_ti(Dpol_b, pol_b, 1);
  Dpol_b[0] &= mask;

  fppol64_t slope; // (a'(omega) - b'(omega)*omega) / b(omega).
  fppol64_t aux;
  fppol64_mul_ti(aux, Dpol_b, 1);
  fppol64_add(aux, aux, Dpol_a);
  fppol64_rem(slope, aux, seven);
  fppol64_rem(aux, pol_b, seven);
  // aux <- 1/aux mod t^2+t+1
  if (aux[0] == (uint64_t)2U)
      aux[0] = (uint64_t)3U;
  else if (aux[0] == (uint64_t)3U)
      aux[0] = (uint64_t)2U;
  fppol64_mul(slope, slope, aux);
  fppol64_rem(slope, slope, seven);

  // Three cases:
  //   slope = 0    ->   [ e-1, 1 ]
  //   slope = 1    ->   [ 1, e-1 ]
  //   slope = other ->  [ 1, 1 ]   and in that case, e = 2.

  if (slope[0] == (uint64_t)0U)
  {
    exp_above[0] = e-1;
    exp_above[1] = 1;
  }
  else if (slope[0] == (uint64_t)1U)
  {
    exp_above[0] = 1;
    exp_above[1] = e-1;
  }
  else
  {
    exp_above[0] = 1;
    exp_above[1] = 1;
    ASSERT_ALWAYS (e == 2);
  }



#else
  #ifndef FOR_FFS
  fprintf (stderr, "Error, the relation with a=%"PRId64" b=%"PRIu64" contains a"
                   " bad ideals (above p=%lu)\n", a, b, p);
  #else
  fprintf (stderr, "Error, the relation with a=%"PRIx64" b=%"PRIx64" contains a"
                   " bad ideals (above p=%lu)\n", (uint64_t) a, b, p);
  #endif
  fprintf (stderr, "but the function handling the bad ideals does not exist.\n"
                   "See filter/filter_badideals.h\n");
  ASSERT_ALWAYS(0);
#endif
}
