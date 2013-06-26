/* Do not check that p is really a bad ideal, must call renumber_is_bad before */
/* The renumbering index should be written in s */
inline void
handle_bad_ideals (prime_t *above, int64_t a, uint64_t b, unsigned long p, int e)
{
/****for a p59 ******/
  if (p == 2) 
  {
    // 0 and 1 are the two ideals above p=2 r=0
    above[0].h = 0;
    above[1].h = 1;

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
      above[0].e = 1;
      above[1].e = 2;
    }
    else if (v == 2) 
    {
      ASSERT_ALWAYS(e == 5);
      above[0].e = 2;
      above[1].e = 3;
    }
    else if (v == 6) 
    {
      ASSERT_ALWAYS(e >= 6);
      above[0].e = e-3;
      above[1].e = 3;
    }
    else
      ASSERT_ALWAYS(0);
  }
  if (p == 3) 
  {
    // 2 and 3 are the two ideals above p=3 r=0
    above[0].h = 2;
    above[1].h = 3;

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
      above[0].e = 1;
      above[1].e = 1;
    }
    else if (v == 3) 
    {
      ASSERT_ALWAYS(e >= 2);
      above[0].e = 1;
      above[1].e = e-1;
    }
    else if (v == 6) 
    {
      ASSERT_ALWAYS(e >= 2);
      above[0].e = e-1;
      above[1].e = 1;
    }
    else
      ASSERT_ALWAYS(0);
  }
  else
    ASSERT_ALWAYS(0);
}


#if 0
  /******For a c155 ***********/
  if (p == 2) 
  {
    // 0 and 1 are the two ideals above p=2 r=0
    above[0].h = 0;
    above[1].h = 1;
    // 1 has always exponent 1
    above[1].e = 1;
    
    //hack: this compute a/b % 16 (not very efficient)
    unsigned long r = findroot (a, b, 16);
    if (r == 4 || r == 6 || r == 10 || r == 14)
    {
      ASSERT_ALWAYS(e == 2); //Is it correct?
      above[0].e = 1;
    }
    else if (r == 0 || r == 8)
    {
      ASSERT_ALWAYS(e == 3); //Is it correct?
      above[0].e = 2;
    }
    else if (r == 4)
    {
      ASSERT_ALWAYS(e == 4); //Is it correct?
      above[0].e = 3;
    }
    else if (r == 12)
    {
      above[0].e = e-1;
    }
    else
      ASSERT_ALWAYS(0);
  }
  else if (p == 5) 
  {
    // 2 and 3 are the two ideals above p=5 r=2
    above[0].h = 2;
    above[1].h = 3;

    //hack: this compute a/b % 16 (not very efficient)
    unsigned long r = findroot (a, b, 25);
    if (r == 5 || r == 7 || r == 22)
    {
      above[0].e = 1;
      above[1].e = 1;
    }
    else if (r == 12)
    {
      above[0].e = 1;
      above[1].e = e-1;
    }
    else if (r == 17)
    {
      above[0].e = e-1;
      above[1].e = 1;
    }
    else
      ASSERT_ALWAYS(0);
  }
  else
    ASSERT_ALWAYS(0);
#endif
