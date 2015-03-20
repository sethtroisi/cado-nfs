#!/usr/bin/env python3
#
# Example for the p59 in tests/test_full_p59
#
#   ./scripts/descent.py		\
#   --target 15384226356769205532939362866999574621778419520781718884300  \
#   --db ~/Local/p59/p59.db	        \
#   --cadobindir $PWD/build/`hostname`	\
#   --init-I 9		                \
#   --init-ncurves 20	                \
#   --init-lpb 24	                \
#   --init-lim 1000000	                \
#   --init-mfb 40	                \
#   --init-tkewness 8388608		\
#   --I 11		\
#   --lpb0 22		\
#   --lpb1 22		\
#   --lim0 50000	\
#   --lim1 100000	\
#   --mfb0 40		\
#   --mfb1 40		\
#   --descent-hint ~/Local/p59/p59.hint
#
# where the hintfile must not be empty, it should contain at least one
# line, as in e.g:
#
#          22@0 0.0038 0.9940 I=11 50000,22,24,1.2 100000,22,24,1.2
#          22@1 0.0038 0.9995 I=11 50000,22,24,1.2 100000,22,24,1.2
#          23@0 0.0038 0.9945 I=11 50000,22,24,1.2 100000,23,43,1.5
#          23@1 0.0038 0.9990 I=11 50000,22,24,1.2 100000,22,24,1.2
#          24@0 0.0038 0.9945 I=11 50000,22,24,1.2 100000,24,45,1.5


# TODO: do we want to have default values, as done here, for the --init-*
# arguments ? I would say no.

# TODO: of course this is currently kludgy, and does not properly wield
# the power of the cadofactor python programs. I need to understand that
# stuff better.

# FIXME: we must be resilient to the case where a prime appears twice,
# once rational, once algebraic ; or two algebraic, but with a different
# prime ideal above. This means we must have the extraprimes list richer.

# TODO: this currently fails when the target is ridiculously small, as
# target=1009 for instance -- see why.

# TODO: keep las awake, if needed.

# TODO: make output less verbose.

import os
import sqlite3
import subprocess
import sys
import argparse
import re
import math
import tempfile
import shutil
import functools
import itertools

# This gives the boring info about the file names for everything, and
# the boilerplate arguments to be fed to binaries.
class GeneralClass(object):

    def declare_args(parser):
        parser.add_argument("--no-wipe",
                help="Keep working files",
                action="store_true")
        parser.add_argument("--datadir",
                help="cadofactor working directory",
                type=str)
        parser.add_argument("--prefix",
                help="project prefix",
                type=str)
        parser.add_argument("--db",
                help="SQLite db file name",
                type=str)
        # the next few are optional file names
        parser.add_argument("--tmpdir", help="Temporary working directory")
        parser.add_argument("--cadobindir",
                help="Cado build directory",
                required=True,
                type=str)
#        parser.add_argument("--todofile",
#                help="Output primes to this toto file",
#                type=str)
        parser.add_argument("--poly",
                help="Polynomial file",
                type=str)
        parser.add_argument("--renumber",
                help="Renumber file",
                type=str)
        parser.add_argument("--badidealinfo",
                help="Badideal info file",
                type=str)
        parser.add_argument("--fb1",
                help="Factor base file for the algebraic side",
                type=str)
        parser.add_argument("--log",
                help="File with known logs",
                type=str)
        # This one is a bit special, we'll create it if it is missing.
        parser.add_argument("--debug-renumber",
                help="File given by debug_renumber",
                type=str)
        # This one applies to both las in the initial step, and
        # reconstructlog in the final step
        parser.add_argument("--threads",
                help="Number of threads to use",
                type=int, default=4)
        # the arguments below are really better fetched form the
        # database.
        parser.add_argument("--ell", help="Group order (a.k.a. ell)")
        parser.add_argument("--nsm0", help="Number of SM on side 0")
        parser.add_argument("--nsm1", help="Number of SM on side 1")
        parser.add_argument("--smexp0", help="SM exponent on side 0")
        parser.add_argument("--smexp1", help="SM exponent on side 1")
        # Those are used both for the middle and lower levels of the
        # descent.
        for side in range(2):
            parser.add_argument("--lpb%d" % side,
                    help="Default large prime bound on side %d" % side,
                    required=True,
                    type=int)

    def __init__(self, args):
        self._conn = None
        self.args = args
        if bool(args.db) == bool(args.prefix and args.datadir):
            raise ValueError("Either --db or the combo --prefix + --datadir must be specified")
        if args.tmpdir:
            self._tmpdir = args.tmpdir
            # do mkdir ???
        else:
            self._tmpdir = tempfile.mkdtemp(dir="/tmp")
        self.hello()

    
    def __connect(self):
        if args.db and not self._conn:
            self._conn = sqlite3.connect(args.db)

    def __getdb(self, query):
        if not args.db:
            return None
        self.__connect()
        self._cursor = self._conn.cursor()
        self._cursor.execute(query)
        v=self._cursor.fetchone()
        self._cursor.close()
        return v

    def __getfile(self, shortname, typical, table, key):
        try:
            v=self.args.__dict__[shortname]
            if v:
                return v
        except KeyError:
            pass
        if args.db:
            v=self.__getdb("select value from %s where key='%s'" % (table, key))
            if v is not None and len(v) > 0:
                return os.path.join(os.path.dirname(args.db), v[0])
        elif args.datadir and args.prefix:
            return os.path.join(args.datadir, args.prefix + "." + typical)
        raise ValueError("no %s file known" % shortname)

    def __getarg(self, shortname, table, key):
        try:
            v=self.args.__dict__[shortname]
            if v:
                return v
        except KeyError:
            pass
        if args.db:
            v=self.__getdb("select value from %s where key='%s'" % (table, key))
            if v is not None and len(v) > 0:
                return v[0]
        raise ValueError("no %s parameter known" % shortname)

    def prefix(self):
        if args.prefix:
            return args.prefix
        else:
            return os.path.basename(args.db).split('.')[0]

    def datadir(self):
        if args.datadir:
            return args.datadir
        else:
            return os.path.dirname(args.db)

    def poly(self):
        return self.__getfile("poly", "polyselect2.poly", "polyselect2", "polyfilename")
    def renumber(self):
        return self.__getfile("renumber", "freerel.renumber.gz", "freerel", "renumberfilename")
    def debug_renumber(self):
        if args.debug_renumber:
            return args.debug_renumber
        else:
            return os.path.join(self.datadir(), self.prefix() + ".debug_renumber.txt")

    def log(self):
        return self.__getfile("log", "reconstructlog.dlog", "reconstructlog", "dlog")
    def badidealinfo(self):
        return self.__getfile("badidealinfo", "magmanmbrthry.badidealinfo", "magmanmbrthry", "badinfofile")
    def fb1(self):
        return self.__getfile("fb1", "factorbase.roots.gz", "factorbase", "outputfile")
    def ell(self):
        return int(self.__getarg("ell", "magmanmbrthry", "ell"))
    def nmaps0(self):
        return int(self.__getarg("nmaps0", "magmanmbrthry", "nmaps0"))
    def nmaps1(self):
        return int(self.__getarg("nmaps1", "magmanmbrthry", "nmaps1"))
    def smexp0(self):
        return int(self.__getarg("smexp0", "magmanmbrthry", "smexp0"))
    def smexp1(self):
        return int(self.__getarg("smexp1", "magmanmbrthry", "smexp1"))
    def lpb0(self):
        return args.lpb0
    def lpb1(self):
        return args.lpb1
    def tmpdir(self):
        return self._tmpdir
    def threads(self):
        return int(args.threads)
    def poly_data(self):
        d={}
        with open(self.poly(), "r") as file:
            for line in file:
                if re.match("^\s*#", line):
                    continue
                if re.match("^\s*$", line):
                    continue
                key,value=line.split(":")
                key = key.strip()
                foo = re.match("^([cY])(\d+)$", key)
                if foo:
                    s,i=foo.groups()
                    if s not in d:
                        d[s]=[]
                    while int(i) >= len(d[s]):
                        d[s]+=[None]
                    d[s][int(i)]=value.strip()
                else:
                    d[key] = value.strip()
        return d

    def p(self):
        d=self.poly_data()
        return int(d["n"])

    def rational_poly():
        d=self.poly_data()
        return [int(x) for x in d["Y"]]
    def algebraic_poly():
        d=self.poly_data()
        return [int(x) for x in d["c"]]

    def __del__(self):
        if not self.args.tmpdir and not self.args.no_wipe:
            shutil.rmtree(self.tmpdir())
        if self._conn:
            self._conn.close()

    def las_bin(self):
        return os.path.join(args.cadobindir, "sieve", "las")
    def debug_renumber_bin(self):
        return os.path.join(args.cadobindir, "misc", "debug_renumber")
    def dup2_bin(self):
        return os.path.join(args.cadobindir, "filter", "dup2")
    def reconstructlog_bin(self):
        return os.path.join(args.cadobindir, "filter", "reconstructlog-dl")

    def lasMiddle_base_args(self):
        # TODO add threads once it's fixed.
        s=[
            self.las_bin() + "_descent",
            "--allow-largesq",
            "--never-discard",  # useful for small computations.
            "--renumber", self.renumber(),
            "--log", self.log(),
            "--fb1", self.fb1(),
            "--poly", self.poly(),
          ]
        return [ str(x) for x in s ]
    def reconstructlog_final_base_args(self):
        s=[
            self.reconstructlog_bin(),
            "-gorder", self.ell(),
            "-sm0", self.nmaps0(),
            "-sm1", self.nmaps1(),
            "-smexp0", self.smexp0(),
            "-smexp1", self.smexp1(),
            "-poly", self.poly(),
            "-logformat", "reconstruct",
            "-mt", self.threads()
          ]
        return [ str(x) for x in s ]

    # There's no las_init_base_args, since DescentUpperClass uses only
    # its very own arguments.

    def hello(self):
        print("Working in GF(p), p=%d" % self.p())
        print("Subgroup considered in GF(p)^* has size %d" % self.ell())
        print("prefix is %s" % self.prefix())


# We need this in order to see which are the rational primes we need to
# descend. We do not care about the algebraic primes, here/
class RatLogBase(object):
    def __init__(self, general):
        self.known=set()
        with open(general.log(),'r') as file:
            for line in file:
                foo = re.match("^(\w+) (\w+) (\w+) rat (\d+)$", line)
                if foo:
                    self.known.add(int(foo.groups()[1], 16))
        print("Found %d known rational logs in %s" %(len(self.known), general.log()))
    def has(self, p):
        return p in self.known

class DescentUpperClass(object):
    def declare_args(parser):
        c = " (specific for the descent bootstrap)"
        parser.add_argument("--init-tkewness",
                help="Tkewness"+c,
                type=int,
                default=2**30)
        parser.add_argument("--init-lim",
                help="Factor base bound"+c,
                default=2**26)
        parser.add_argument("--init-lpb",
                help="Large prime bound"+c,
                default=64)
        parser.add_argument("--init-mfb",
                help="Cofactor bound"+c,
                default=100)
        parser.add_argument("--init-ncurves",
                help="ECM effort in cofactorization"+c,
                default=80)
        parser.add_argument("--init-I",
                help="Sieving range"+c,
                default=14)

    def __init__(self, general, known, args):
        self.general = general
        self.known = known

        self.tkewness = int(args.init_tkewness)
        self.lim      = int(args.init_lim)
        self.lpb      = int(args.init_lpb)
        self.mfb      = int(args.init_mfb)
        self.ncurves  = int(args.init_ncurves)
        self.I        = int(args.init_I)

    def __myxgcd(self, a, b, T):
        assert type(a) == int
        assert type(b) == int
        assert type(T) == int
        ainit = a
        binit = b
        bound = math.floor(math.sqrt(b*T))
        x = 0
        lastx = 1
        y = 1
        lasty = 0
        while abs(b) > bound:
            q = a // b
            r = a % b
            a = b
            b = r
            newx = lastx - q*x
            lastx = x
            x = newx
            newy = lasty - q*y
            lasty = y
            y = newy
        return [ [ b, x ], [ a, lastx ] ]

    def do_descent(self, z):
        p = general.p()
        gg = self.__myxgcd(z, p, self.tkewness)
        tmpdir = general.tmpdir()
        prefix = general.prefix() + "-descent_init."

        polyfilename = os.path.join(tmpdir, prefix + "poly")
        with open(polyfilename, 'w') as f:
            f.write("n: %d\n" % p)
            f.write("skew: 1\n")
            f.write("c1: %d\n" % gg[0][0])
            f.write("c0: %d\n" % gg[1][0])
            f.write("Y1: %d\n" % gg[0][1])
            f.write("Y0: %d\n" % gg[1][1])

#       These are not read anyway, and the code in las is much faster.
#        print ("Creating factor bases.")
#        for side in range(2):
#            subprocess.check_call([args.cadobindir + "/sieve/makefb",
#                "-t", str(args.nthreads),
#                "-poly", polyfilename,
#                "-alim", str(args.lim),
#                "-side", str(side),
#                "-out", wdir + str("/desc.fb") + str(side)],
#                stdout=subprocess.DEVNULL)
        
        print ("--- Sieving (initial) ---")
        q0 = self.tkewness
        q1 = self.tkewness + 100000
        # relsfilename = os.path.join(tmpdir, prefix + "rels")

        call_that = [ general.las_bin(),
            "-poly", polyfilename,
            # "-fb0", wdir + str("/desc.fb") + str(0),
            # "-fb1", wdir + str("/desc.fb") + str(1),
            "-lim0", self.lim,
            "-lim1", self.lim,
            "-lpb0", self.lpb,
            "-lpb1", self.lpb,
            "-mfb0", self.mfb,
            "-mfb1", self.mfb,
            "-ncurves0", self.ncurves,
            "-ncurves1", self.ncurves,
            "-I", self.I,
            "-q0", q0,
            "-q1", q1,
            "--exit-early", 2,
            ]
        call_that = [str(x) for x in call_that]
        print("command line:\n" + " ".join(call_that))
        with subprocess.Popen(call_that, stdout=subprocess.PIPE) as las:
            rel = None
            for line in las.stdout:
                line=line.decode("utf-8")
                if line[0] == '#':
                    if (re.match("^# \d+ relation", line)):
                        sys.stdout.write('\n')
                        print(line.rstrip())
                    else:
                        sys.stdout.write('.')
                        sys.stdout.flush()
                    continue;
                else:
                    rel = line.strip()
                    break
        sys.stdout.write('\n')
        if not rel:
            raise NameError("No relation found!")
        print("Taking relation %s\n" % rel)
        rel = rel.split(':')
        a,b = [int(x) for x in rel[0].split(',')]

        Num = a*gg[0][0] + b*gg[1][0]
        Den = a*gg[0][1] + b*gg[1][1]
        assert (z*Den-Num) % p == 0

        factNum = [ int(x, 16) for x in rel[2].split(',') ]
        factDen = [ int(x, 16) for x in rel[1].split(',') ]
        print(Num, Den, factNum, factDen)

        assert(abs(Num) == functools.reduce(lambda x,y:x*y,factNum,1))
        assert(abs(Den) == functools.reduce(lambda x,y:x*y,factDen,1))

        todofilename = os.path.join(tmpdir, prefix + "todo")
        with open(todofilename, "w") as f:
            for q in factNum + factDen:
                if known.has(q):
                    continue
                logq = math.ceil(math.log(q, 2))
                print("Will do further descent for %d-bit rational prime %d"
                        % (logq, q))
                # las can understand when the rational root is missing
                f.write("0 %d\n" % q)
        return todofilename, [Num, Den, factNum, factDen]


class DescentMiddleClass(object):
    def declare_args(parser):
        # TODO: import default values from the sieving parameters.
        parser.add_argument("--descent-hint",
            help="Hintfile for the descent",
            required=True,  # TODO: fall back on default values.
            )
        parser.add_argument("--I",
            help="Default value for I (must match hint file)",
            required=True,
            type=int)
        for side in range(2):
            parser.add_argument("--mfb%d" % side,
                    help="Default cofactor bound on side %d" % side,
                    required=True,
                    type=int)
            parser.add_argument("--lim%d" % side,
                    help="Default factor base bound on side %d (must match hint file)" % side,
                    required=True,
                    type=int)

    def __init__(self, general, args):
        self.general = general
        self.args = args
        # We need to do some safety checking
        values_I=set()
        values_lim0=set()
        values_lim1=set()
        values_I.add(args.I)
        values_lim0.add(args.lim0)
        values_lim1.add(args.lim1)
        with open(args.descent_hint, 'r') as file:
            for line in file:
                if re.match("^\s*#", line):
                    continue
                if re.match("^\s*$", line):
                    continue
                line = line.strip()
                foo = re.match("^.*I=(\d+)\s+(\d+),[\d.,]+\s+(\d+),[\d.,]+$",
                        line)
                if not foo:
                    print("Warning, parse error in hint file at line:\n"
                            + line);
                    continue
                I,lim0,lim1 = foo.groups()
                values_I.add(int(I))
                values_lim0.add(int(lim0))
                values_lim1.add(int(lim1))
        if len(values_lim0)>1:
            raise ValueError("lim0 values should match between cmdline and hint file")
        if len(values_lim1)>1:
            raise ValueError("lim1 values should match between cmdline and hint file")
        if len(values_I)>1:
            raise ValueError("I values should match between cmdline and hint file")
        print("Consistency check for las_descent passed");
        print("\tI=%d" % values_I.pop())
        print("\tlim0=%d" % values_lim0.pop())
        print("\tlim1=%d" % values_lim1.pop())


    def do_descent(self, todofile):
        tmpdir = general.tmpdir()
        prefix = general.prefix() + "-descent_middle."
        f = open(todofile, 'r')
        ntodo = len(list(f))
        f.close()
        print ("--- Sieving (middle, %d rational primes) ---" % ntodo)
        s=general.lasMiddle_base_args()
        if args.descent_hint:
            s += [ "--descent-hint-table", args.descent_hint ]
        s += [
                "--I", self.args.I,
                "--lim0", self.args.lim0,
                "--lpb0", general.lpb0(),
                "--mfb0", self.args.mfb0,
                "--lim1", self.args.lim1,
                "--lpb1", general.lpb1(),
                "--mfb1", self.args.mfb1,
             ]
        s += [ "--todo", todofile ]
        relsfilename = os.path.join(tmpdir, prefix + "rels")
        s += [ "--out", relsfilename ]
        call_that=[str(x) for x in s]
        print("command line:\n" + " ".join(call_that))
        subprocess.check_call(call_that)
        return relsfilename

class DescentLowerClass(object):
    def declare_args(parser):
        pass
    def __init__(self, general, args):
        self.general = general
        self.args = args
    def _last_renumber_index(self):
        # In the renumber file, the line number is not exactly the index.
        # Let us rely on debug_renumber's output to find the last used index
        # in the renumber table.
        if not os.path.exists(general.debug_renumber()):
            print("--- generating %s ---"% general.debug_renumber())
            call_that = [ general.debug_renumber_bin(),
                            "-poly", general.poly(),
                            "-renumber", general.renumber(),
                        ]
            call_that = [str(x) for x in call_that]
            outfile = open(general.debug_renumber(), 'w')
            subprocess.check_call(call_that, stdout=outfile)
            outfile.close()

        last_renumber_line = None

        with subprocess.Popen(["tail", "-10", general.debug_renumber()],
                stdout=subprocess.PIPE) as p:
            for line in p.stdout:
                ll = line.decode("utf-8")
                if ll[0] == 'i':
                    last_renumber_line = ll.rstrip()
        return int(last_renumber_line.split()[0].split('=')[1], 16)

    def do_descent(self, relsfile, initial_split):
        args = parser.parse_args()
        tmpdir = general.tmpdir()
        prefix = general.prefix() + "-descent_lower."

        # Read descent relations
        noncomment=lambda x: not re.match("^#",x)
        with open(relsfile, 'r') as file:
            descrels = list(filter(noncomment, file))
        nrels = len(descrels)
        print ("--- Final reconstruction (from %d relations) ---" % nrels)

        lpb = [ general.lpb0(), general.lpb1() ]

        # Create fake relations: remove large primes, because those are
        # not in the renumber table.
        #
        # FIXME -- we must create proper renumber table entries for extra
        # primes, otherwise the birthday paradox will kill us.
        fakerels = []
        extraprimes = []
        for rel in descrels:
            r = rel.split(':')
            ss = [ [], [] ]
            extra = [ [], [] ]
            for side in range(2):
                for p in r[side+1].strip().split(','):
                    if int(p, 16) >> lpb[side]:
                        extra[side].append(p)
                    else:
                        ss[side].append(p)
            fake = ":".join([r[0]] + [",".join(x) for x in ss])
            fakerels.append(fake)
            extraprimes.append(extra)
        setextraprimes = set(sum(sum(extraprimes,[]),[]))
        listextraprimes = list(setextraprimes)
        numlist = [ int(p, 16) for p in listextraprimes ]
        numlist.sort()
        listextraprimes = [ hex(p)[2:] for p in numlist ]
        print ("They include %s primes above lpb:" % len(listextraprimes),
                ", ".join(listextraprimes))
        fakefilename = os.path.join(tmpdir, prefix + "fakerels")
        with open(fakefilename, 'w') as f:
            for rel in fakerels:
                f.write(rel + '\n')

        # Call dup2 to renumber those relations to get them renumbered
        # It works in place, so we do a copy first.
        fakefilename2 = os.path.join(tmpdir, prefix + "fakerels-renumbered.txt")

        shutil.copyfile(fakefilename, fakefilename2)

        print("--- renumbering descent relations ---")
        call_that = [ general.dup2_bin(),
                        "-dl",
                        "-nrels", nrels,
                        "-renumber", general.renumber(),
                        "-badidealinfo", general.badidealinfo(),
                        fakefilename2
                    ]
        call_that = [str(x) for x in call_that]
        subprocess.check_call(call_that, stderr=subprocess.DEVNULL)

        last_index = self._last_renumber_index()
        print("--- last renumber index is %d ---" % last_index)



        # Now, we can assign fresh indices to new primes
        # We fix the relations and the create a new renumber table
        # accordingly.
        dictextraprimes = {}
        for i in range(len(listextraprimes)):
            dictextraprimes[listextraprimes[i]] = last_index+1+i
        print(dictextraprimes)

        descentrelsfilename = os.path.join(tmpdir, prefix + "descentrels-renumbered.txt")
        # relations
        with open(descentrelsfilename, 'w') as outfile:
            # add first line like in purged.gz, just to give the number
            # of relations to reconstructlog.
            outfile.write("# %d 0 0\n" % nrels)
            extra = iter(extraprimes)
            with open(fakefilename2, 'r') as infile:
                for rel in infile:
                    rel = rel.rstrip()
                    for p in sum(next(extra),[]):
                        rel += "," + hex(dictextraprimes[p])[2:]
                    outfile.write(rel + "\n")

        newrenumberfilename = os.path.join(tmpdir, prefix + "new-renumber.txt")
        print("--- generating temporary %s ---" % newrenumberfilename)
        # we now have to create a new temp file which contains the
        # renumber table + our extra crap appended.
        newlpb = math.ceil(math.log(int(listextraprimes[-1], 16), 2))
        with open(newrenumberfilename, 'w') as file:
            # decompress the old one.
            with subprocess.Popen([ "gzip", "-dc", general.renumber() ],
                    stdout=subprocess.PIPE) as old:
                # tweak the header
                header = old.stdout.readline().decode("utf-8").strip().split()
                header[-1] = str(newlpb)
                header[-2] = str(newlpb)
                file.write(" ".join(header) + "\n")
                # leave the bulk of the file unchanged
                for x in old.stdout:
                    file.write(x.decode("utf-8"))
                # add our crap
                for p in listextraprimes:
                    file.write(hex(int(p, 16)+1)[2:]+'\n')

        newdlogfilename = os.path.join(tmpdir, prefix + "new-dlog.txt")
        print("--- generating temporary %s ---" % newdlogfilename)
        # Prepare the call to reconstructlog
        # Need to modifiy the last lines of known logs that contain the
        # virtual logs of the SMs
        nsm = general.nmaps0() + general.nmaps1()
        with open(newdlogfilename, 'w') as outfile:
            count_sm=0
            with open(general.log(), 'r') as infile:
                for line in infile:
                    if "SM" in line:
                        count_sm += 1
                        ll = line.split()
                        assert len(ll) == 5
                        ll[0] = hex(int(ll[0], 16) + len(listextraprimes))[2:]
                        line = " ".join(ll) + "\n"
                    outfile.write(line)
            assert nsm == count_sm

        # create an empty relsdel file
        fakerelsdelfilename = os.path.join(tmpdir, prefix + "fake-relsdel.txt")
        open(fakerelsdelfilename,'w').close()

        resultfilename = os.path.join(tmpdir, prefix + "result-dlog.txt")
        s = general.reconstructlog_final_base_args()
        s+= [
            "-nrels", nrels,
            "-renumber", newrenumberfilename,
            "-log", newdlogfilename,
            "-purged", descentrelsfilename,
            "-relsdel", fakerelsdelfilename,
            "-out", resultfilename,
            ]
        call_that=[str(x) for x in s]
        print("command line:\n" + " ".join(call_that))
        subprocess.check_call(call_that)

        ## Deduce the log of the target, at last...
        Num, Den, factNum, factDen = initial_split
        # need to grep these rational primes in our log table. We do not
        # care about the index, really, but match with the rest.
        wanted={ hex(int(p))[2:]:None for p in factNum+factDen }
        # oh, and by the way, let's forcibly add 2 and 3 to the list of
        # wanted logs. This way, the user will be able to check the
        # consistency easily
        wanted["2"]=None
        wanted["3"]=None
        with open(resultfilename, 'r') as file:
            for line in file:
                foo=re.match("^\w+ (\w+) 0 rat (\d+)$", line)
                if foo:
                    pp = foo.groups()[0]
                    if pp in wanted:
                        wanted[pp]=int(foo.groups()[1])
        log_target = 0
        for p in factNum:
            pp = hex(int(p))[2:]
            log_target = log_target + wanted[pp]
        for p in factDen:
            pp = hex(int(p))[2:]
            log_target = log_target - wanted[pp]
        p=general.p()
        ell=general.ell()
        print("# p=%d" % p)
        print("# ell=%d" % ell)
        print("log(2)=%d" % wanted["2"])
        print("log(3)=%d" % wanted["3"])
        print("# target=%s" % args.target)
        print("log(target)=%d" % log_target)
        print("# check line: " +
                "((GF(p)!2)^%d/(GF(p)!%s)^%d)^((p-1)div %d) where p is %d;"
                % (log_target, args.target, wanted["2"], ell, p))


if __name__ == '__main__':

    print("THIS IS AN EXPERIMENTAL SCRIPT ONLY ! STILL AT WORK.")

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Descent initialization for DLP")

    # Required
    parser.add_argument("--target", help="Element whose DL is wanted",
            type=int, required=True)

    GeneralClass.declare_args(parser)
    DescentUpperClass.declare_args(parser)
    DescentMiddleClass.declare_args(parser)
    DescentLowerClass.declare_args(parser)

    args = parser.parse_args()

    general = GeneralClass(args)
    known = RatLogBase(general)
    init = DescentUpperClass(general, known, args)
    middle = DescentMiddleClass(general, args)
    lower = DescentLowerClass(general, args)

    todofile, initial_split = init.do_descent(int(args.target))
    relsfile = middle.do_descent(todofile)
    lower.do_descent(relsfile, initial_split)

