#!/usr/bin/env python3

import sys
import glob
import re
import subprocess

class GnuplotFile(object):
    def __init__(self, outputfile, size=(1024, 768), yrange=(55, 59)):
        self.lines = [
            "set terminal png size %d,%d enhanced font 'Helvetica,20'\n" % size,
            "set output '%s'\n" % outputfile,
            "set yrange [%d:%d]\n" % yrange
        ]
        self.is_first = True

    def __str__(self):
        return "".join(self.lines) + "\n"

    def append_plot(self, filename, title):
        prefix = "plot" if self.is_first else ","
        self.lines.append("%s '%s' title '%s'" % (prefix, filename, title))
        self.is_first = False

def grep_files(path, regex=None):
    if not regex is None:
        compiled_re = re.compile(regex)
    for filename in glob.glob(path):
        with open(filename, "r") as input_file:
            for line in input_file:
                if regex is None or compiled_re.search(line):
                    yield line

def make_polyselect_filename(directory, adstart, adlen):
    """ Returns the polyselect2l output file name for a given admin -- admax range """
    return '%s/c200.polyselect1.*.%s-%s' % (directory, adstart, adstart + adlen)

def make_polyselect_filenames(directory, admin, admax, adstep, adlen):
    return ((adstart, make_polyselect_filename(directory, adstart, adlen))
            for adstart in range(admin, admax, adstep))

def grep_poly_files(directory, admin, admax, adstep, adlen, regex=None):
    for (adstart, path) in make_polyselect_filenames(directory, admin, admax, adstep, adlen):
        for line in grep_files(path, regex):
            yield (adstart, line)

def get_nr_poly(directory, admin, admax, adstep, adlen, max_nr_poly=0):
    total = 0
    max_ad = None
    for (adstart, line) in grep_poly_files(directory, admin, admax, adstep, adlen, regex="^# Stat: tried"):
        # The max_ad value that gets returned means that we can use ad in [0, max_ad] to get no more than
        # max_nr_poly polynomials. We need to update max_ad *before* updating total, because the currently
        # processed file is for the [adstart, adstart + adlen] interval.
        if not max_nr_poly is None and total < max_nr_poly:
            max_ad = adstart
        total += int(line.split()[6])
    return (total, max_ad)

def get_time(directory, admin, admax, adstep, adlen):
    total = 0.
    for (adstart, line) in grep_poly_files(directory, admin, admax, adstep, adlen, regex="^# Stat: total phase took"):
        total += float(line.split()[5].strip("s"))
    return total

def make_best_lognorms(P, nq, nrkeep, admin, admax, adstep, adlen,
                       gnuplot_file=None):
    directory = "P_%d.nq_%d" % (P, nq)
    output_filename = "%s/lognorms.%d-%d" % (directory, admin, admax)
    lognorms = []
    print("Generating list of best %d lognorms in %s from ad=%d up to %d in "
          "file %s" %
          (nrkeep, directory, admin, admax, output_filename))
    for (adstart, line) in grep_poly_files(directory, admin, admax, adstep, adlen, regex="^# lognorm"):
        lognorms.append(float(line.split()[2].strip(",")))
    lognorms.sort()
    lognorms = lognorms[0:nrkeep]
    print("Best lognorm was %f, %d-th best was %f" % (lognorms[0], len(lognorms), lognorms[-1]))
    with open(output_filename, "w") as output_file:
        output_file.write("\n".join(map(str, lognorms)))
    if gnuplot_file:
        gnuplot_file.append_plot(output_filename, title="P=%d, ad in [%d, %d]"
                                 % (P, admin, admax))

def run_gnuplot(plot):
    filename = "input.gnuplot"
    with open(filename, "w") as output_file:
        output_file.write(str(plot))
    subprocess.check_call(["gnuplot", filename])

def process_one_nq(nq, nrkeep, admin, admax, adstep, adlen, max_nr_poly):
    plot_allruns = GnuplotFile("graph.best_%d.nq_%d.allruns.png" % (nrkeep, nq))
    plot_only13k = GnuplotFile("graph.best_%d.nq_%d.%d.png" % (nrkeep, nq, max_nr_poly))
    for P in P_VALUES:
        directory = "P_%d.nq_%d" % (P, nq)
        (nr_poly, admax_13k) = get_nr_poly(directory, admin, admax, adstep, adlen, max_nr_poly=max_nr_poly)
        for (admax_now, plot_now) in zip([admax, admax_13k],
                                         [plot_allruns, plot_only13k]):
            (nr_poly, _) = get_nr_poly(directory, admin, admax_now, adstep, adlen)
            time = get_time(directory, admin, admax_now, adstep, adlen)
            print("Number of polynomials and time spent in %s up to ad=%d: %d / %g" %
                  (directory, admax_now, nr_poly, time))
            make_best_lognorms(P, nq, nrkeep, admin, admax_now, adstep,
                               adlen, plot_now)

    print("Generating plots")
    run_gnuplot(plot_allruns)
    run_gnuplot(plot_only13k)

assert len(sys.argv) == 7

ADMIN = int(sys.argv[1])
ADMAX = int(sys.argv[2])
ADSTEP = int(sys.argv[3])
ADLEN = int(sys.argv[4])
NRKEEP = int(sys.argv[5])
MAX_NR_POLY = int(sys.argv[6])

NQ_VALUES = [216]
P_VALUES = [10000, 100000, 1000000, 10000000]

for NQ in NQ_VALUES:
    process_one_nq(NQ, NRKEEP, ADMIN, ADMAX, ADSTEP, ADLEN, MAX_NR_POLY)
