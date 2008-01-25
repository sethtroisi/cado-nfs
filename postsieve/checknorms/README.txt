README.txt file updated on Fri Jan 25 2008 by JM

1. What is it?
--------------

  This is a slightly revamped variant of Paul's checknorms program. It has
  been modified to use TIFA instead of ECM. Usage is mostly unchanged:

      checknorms -poly <file.poly> <file.raw.rels> > <file.rels>

  However new options are avaible (type ./checknorms -h).
    
     Usage:
     ------
       ./checknorms [-h] [-v] [-cmult]
                    [-maxnlp <num>] [-mfbr <num>] [-mfba <num>] [-t <num>]
                    -poly <file> <relfile_1> [<relfile_2> ... <relfile_n>]
    
     Mandatory arguments:
     --------------------
       -poly FILE
           CADO polynomial file.
       <relfile_i> FILES
           Space-separated list of relation files obtained from CADO siever.
    
     Options:
     --------
       -v
           Turn verbose mode on.
       -h
           Print this help.
       -cmult
           Take into account multiplicities of residues' factors.
       -maxnlp NUM
           Maximum number of large primes to allow.
           Multiplicities are taken into account if the option -cmult is given.
           Default: 8
       -mfbr NUM
           Bound for rational residues (in bits).
           Default: use value from polynomial file
       -mfba NUM
           Bound for algebraic residues (in bits).
           Default: use value from polynomial file
       -t NUM
           Bound for largest prime for trial division.
           Default: 100
    
  WARNING: Please note that the '-factorall' option has been removed!
  
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
