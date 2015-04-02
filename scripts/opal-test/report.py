#!/usr/bin/env python3

import sys
import re
from math import log, sqrt
try:
  from collections.abc import Iterable
except ImportError:
  from collections import Iterable

RE_FP = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
CAP_FP = "(%s)" % RE_FP
REGEXES = {"cap_fp" : CAP_FP}

PATTERN_SIEVE_SQ = re.compile(r"# Sieving algebraic q=(\d+);")
PATTERN_SQ = re.compile(r"# Average J=(\d+) for (\d+) special-q's, max bucket fill {cap_fp}".format(**REGEXES))
PATTERN_CPUTIME = re.compile(r"# Total cpu time {cap_fp}s .norm {cap_fp}\+{cap_fp}, sieving {cap_fp} .{cap_fp} \+ {cap_fp} \+ {cap_fp}., factor {cap_fp}.".format(**REGEXES))
PATTERN_ELAPSEDTIME = re.compile("# Total elapsed time {cap_fp}s, per special-q {cap_fp}s, per relation {cap_fp}s".format(**REGEXES))
PATTERN_REPORTS = re.compile(r"# Total (\d+) reports .{cap_fp}s/r, {cap_fp}r/sq.".format(**REGEXES))
PATTERN_DUPE = re.compile("# DUPE ")

class NumInt(object):
    """
    >>> n = NumInt()
    >>> n.add(0,1)
    >>> n.add(1,2)
    >>> n.get_value()
    1.5
    >>> n.add(2,3)
    >>> n.get_value()
    4.0
    >>> x = n.interpolate_for_value(3.); 1.6457 < x < 1.6458
    True
    >>> 2.9999 < n.interpolate_at_coord(x) < 3.0001
    True
    """
    def __init__(self):
        # At index 0, the most recent value,coordinate pair.
        # At index 1, the value,coordinate pair before that.
        self.lastvalue = [None, None]
        self.lastcoord = [None, None]
        self.sum = 0
    def trapez_area(self):
        """ Return area of the trapezoidal defined by the last two pairs of coordinates """
        return (self.lastcoord[0] - self.lastcoord[1])*(self.lastvalue[0] + self.lastvalue[1]) / 2.
    def add(self, coord, value):
        self.lastvalue[1], self.lastcoord[1] = self.lastvalue[0], self.lastcoord[0]
        self.lastvalue[0], self.lastcoord[0] = value, coord
        if not self.lastcoord[1] is None:
            assert coord > self.lastcoord[1]
            self.sum += self.trapez_area()
    def get_value(self):
        return self.sum
    def interpolate_for_value(self, value):
        """ Find a coordinate c, greater than the second-last one, such that
        cutting or extending the last trapezoidal up to abscissa c results in
        sum = "value"
        """
        prev_sum = self.sum - self.trapez_area()
        diff = value - prev_sum
        t = (self.lastvalue[0] - self.lastvalue[1]) / (self.lastcoord[0] - self.lastcoord[1])
        # Choose offset x such that, with c0 = c1 + x and v0 = v1 + x * t,
        # prev_sum + (c0 - c1)*(v0 + v1) / 2. = value
        # (c0 - c1)*(v0 + v1) / 2. = diff
        # c0 - c1 = x, v0 + v1 = 2*v1 + x*t
        # x*(2*v1 + x*t) / 2. = diff
        # t/2*x^2 + v1*x - diff = 0
        # x = (-v1 +- sqrt(v1^2 + 2*t*diff)) / t
        # We need only the positive solution
        v1 = self.lastvalue[1]
        disc = v1**2 + 2*t*diff
        if disc < 0:
            sys.stderr.write("discriminant = %f < 0! t = %f, diff = %d, lv=%s, lc=%s\n" % (disc, t, diff, self.lastvalue, self.lastcoord))
        x = (-v1 + sqrt(disc)) / t
        return self.lastcoord[1] + x
    def interpolate_at_coord(self, coord):
        """ Return the sum that would result if the last trapezoidal had been
        cut or extended to abscissa "coord".
        """
        # assert self.lastcoord[1] <= coord <= self.lastcoord[0]
        x = coord - self.lastcoord[1]
        prev_sum = self.sum - self.trapez_area()
        t = (self.lastvalue[0] - self.lastvalue[1]) / (self.lastcoord[0] - self.lastcoord[1])
        v0 = self.lastvalue[1] + x * t
        return prev_sum + x*(self.lastvalue[1] + x * t / 2.)
        

class ListArith(list):
    """
    >>> a = ListArith([1,2,3])
    >>> b = ListArith([3,4,5])
    >>> a + 1
    [2, 3, 4]
    >>> a - 1
    [0, 1, 2]
    >>> a * 2
    [2, 4, 6]
    >>> a + b
    [4, 6, 8]
    >>> b - a
    [2, 2, 2]
    >>> a * b
    [3, 8, 15]
    """
    def __add__(self, other):
        if isinstance(other, Iterable):
            return ListArith([a + b for a,b in zip(self, other)])
        else:
            return ListArith([a + other for a in self])

    def __sub__(self, other):
        if isinstance(other, Iterable):
            return ListArith([a - b for a,b in zip(self, other)])
        else:
            return ListArith([a - other for a in self])

    def __mul__(self, other):
        if isinstance(other, Iterable):
            return ListArith([a * b for a,b in zip(self, other)])
        else:
            return ListArith([a * other for a in self])

    def to_str(self):
        formats = ["%s"] * len(self)
        pat = " ".join(formats)
        return pat % tuple(self)


class LasStats(object):
    def __init__(self):
        self.dupes = 0
        self.J_sum = 0.
        self.nr_sq = 0
        self.max_fill = 0.
        self.cputimes = ListArith([0.] * 8)
        self.eltimes = ListArith([0.] * 3)
        self.reports = 0
        self.relations_int = NumInt()
        self.dupes_int = NumInt()
        self.elapsed_int = NumInt()

    def parse_one_input(self, lines, verbose=False):
        nr_sq = None
        cputimes = None
        eltimes = None
        reports = None
        new_dupes = 0
        first_sq = None
        for line in lines:
            if PATTERN_DUPE.match(line):
                new_dupes += 1
            match = PATTERN_SIEVE_SQ.match(line)
            if match:
                last_sq = int(match.group(1))
                if first_sq is None:
                    first_sq = last_sq
                else:
                    assert first_sq <= last_sq
            match = PATTERN_SQ.match(line)
            if match:
                avg_J = float(match.group(1))
                nr_sq = int(match.group(2))
                new_max_fill = float(match.group(3))
            match = PATTERN_CPUTIME.match(line)
            if match:
                cputimes = list(map(float, match.groups()))
            match = PATTERN_ELAPSEDTIME.match(line)
            if match:
                eltimes = list(map(float, match.groups()))
            match = PATTERN_REPORTS.match(line)
            if match:
                reports = int(match.group(1))
        if nr_sq is None:
            sys.stderr.write("Did not receive value for nr_sq\n")
            return False
        if cputimes is None:
            sys.stderr.write("Did not receive value for cputimes\n")
            return False
        if eltimes is None:
            sys.stderr.write("Did not receive value for eltimes\n")
            return False
        if reports is None:
            sys.stderr.write("Did not receive value for reports\n")
            return False
        self.nr_sq += nr_sq
        self.dupes += new_dupes
        self.J_sum += avg_J * nr_sq
        self.reports += reports
        self.max_fill = max(self.max_fill, new_max_fill)
        self.cputimes += cputimes
        self.eltimes += eltimes
        cputimes_str = self.cputimes.to_str()
        eltimes_str = self.eltimes.to_str()
        sq = (last_sq + first_sq) / 2
        sq_correction = 1./nr_sq/log(sq)
        self.relations_int.add(sq, reports * sq_correction)
        self.dupes_int.add(sq, new_dupes * sq_correction)
        self.elapsed_int.add(sq, eltimes[0]  * sq_correction)
        if verbose:
            names = ("sq", "avgJ", "nr_sq", "sq_sum", "max_fill", "cputimes_str", "elapsed", "elapsed/sq", "elapsed/rel", "reports", "reports/nr_sq", "reports/sqrange", "dupes")
            values = (sq, self.J_sum / self.nr_sq, nr_sq, self.nr_sq, self.max_fill, cputimes_str, eltimes_str, reports, reports/nr_sq, reports * sq_correction, self.dupes)
            print(", ".join( (":".join(map(str, x)) for x in zip(names, values)) ))
        return True

    def parse_one_file(self, filename, verbose=False):
        with open(filename, "r") as f:
            return self.parse_one_input(f, verbose)
     
    def get_rels(self):
        return self.relations_int.get_value()

    def get_dupes(self):
        return self.dupes_int.get_value()

    def get_qmax(self, nr_relations):
            return self.relations_int.interpolate_for_value(nr_relations)

    def get_time(self, nr_relations=None):
        if nr_relations is None:
            return self.elapsed_int.get_value()
        else:
            qmax = self.get_qmax(nr_relations)
            return self.elapsed_int.interpolate_at_coord(qmax)

    def print_stats(self):
        print("Estimated total relations: %f" % self.get_rels())
        print("Estimated total dupes: %f" % self.get_dupes())
        print("Estimated total elapsed time: %f" % self.get_time())
        print("Estimated relations/second: %f" % (self.get_rels() / self.get_time()))


def run():
    stats = LasStats()
    for filename in sys.argv[1:]:
        print("Parsing file %s" % filename)
        stats.parse_one_file(filename, verbose=True)
    stats.print_stats()

if __name__ == '__main__':
    run()
