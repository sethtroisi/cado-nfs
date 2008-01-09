README.txt file updated on Tue Jan 8 2008 by JM

1. What is it?
--------------

  This is a slightly revamped variant of Paul's checknorms program. It has
  been modified to use TIFA instead of ECM. Usage is unchanged:

      checknorms -poly <file.poly> <file.raw.rels> > <file.rels>

  Note that no warning is issued if some prime factor of the residue is larger
  than the factor base bound contrary to Paul's version.

  Also added a '-factorall' option changing the behaviour of the program when
  factoring the norms. What the '-factorall' option does is:
      - Set the tdmax parameter to min(10, tdmax)
      - Attempt to completely factor the norms and accept relations if norm can
        be written as pa^ma * pb^mb where pa and pb are two (possibly
        identical) large primes.

  Note that using the '-factorall' option will slow down the program
  noticeably. My feeling is that _usually_ you don't need to care about it...

2. How to compile it?
---------------------

  To compile it, you will obviously need to compile TIFA (warning: revision
  504 or higher needed!). The provided Makefile has a switch called
  USE_TIFA_FROM_SVN. If set to one (1), the needed TIFA files (library and 
  headers) will be searched in the SVN repository. Note that you will probably 
  have to modify the TIFALIB variable since TIFA's building directory may 
  depends on the host and OS used (this depends on the options chosen when 
  compiling TIFA). If USE_TIFA_FROM_SVN is different from one (1) the needed 
  files will be searched in the directory where TIFA was installed.

  Setting the preprocessor symbol DEBUG to 1 (line 30 of checknorms.c) will
  change the behavior of the program to prompt for user input before reading a
  relation and will print a lots of messages.

3. Not so fast!
---------------

  Finally, a shameless word of warning:

    Don't forget to set the VERBOSE and TIMINGS options to 'no' in TIFA's
    BuildOptions.py file or you'll end up with traces from SQUFOF or SIQS
    in your relation files (and with gigabytes on your hard drive!).
