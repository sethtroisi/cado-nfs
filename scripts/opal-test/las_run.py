from opal.core.io import *
import sys
import subprocess
import re
from math import log
import os

from report import LasStats

def update_existing(a, b):
    """ Update those keys in "a" which also exist in "b" with the values
    from "b"
    """
    for key in set(a) & set(b):
        a[key] = b[key]

# return prime_pi(2^x)
def primepi(x):
    # precomputed values up to x=40 (obtained with Sage)
    l = [0, 1, 2, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512, 6542, 12251, 23000, 43390, 82025, 155611, 295947, 564163, 1077871, 2063689, 3957809, 7603553, 14630843, 28192750, 54400028, 105097565, 203280221, 393615806, 762939111, 1480206279, 2874398515, 5586502348, 10866266172, 21151907950, 41203088796]
    return l[x]

# return the best time so far (or 1e308 at the beginning)
def get_best_time(file):
    try:
        f = open(file, "r")
        best_time = float(f.read())
        f.close()
    except (OSError, IOError):
        best_time = 1e308
    return best_time

def set_best_time(file, t):
    f = open(file, "w")
    f.write(str(t) + "\n")
    f.close()

def run(param_file, problem):
    "Run las with given parameters until the required number of relations is found."

    build_dir = os.environ["CADO_BUILD"]

    makefb = "%s/sieve/makefb" % build_dir
    las = "%s/sieve/las" % build_dir

    las_params = {
        "I": 11,
        "poly": "c59.polyselect2.poly",
        "fb": "c59.factorbase.roots.gz",
        "rlim": 50000,
        "alim": 100000,
        "lpbr": 22,
        "lpba": 22,
        "mfbr": 22,
        "mfba": 22,
        "rlambda": 1.2,
        "alambda": 1.2,
        "ncurves0": 6,
        "ncurves1": 6,
        "t": 1 # number of threads for las
    }
    makefb_params = {
        "poly" : las_params["poly"],
        "alim": las_params["alim"],
        "maxbits": las_params["I"] - 1
    }

    params = read_params_from_file(param_file)

    # Update parameters for las
    update_existing(las_params, params)
    las_params["mfbr"] = max(las_params["mfbr"], las_params["lpbr"])
    las_params["mfba"] = max(las_params["mfba"], las_params["lpba"])
    
    las_params["alambda"] = 1.0 * las_params["mfba"] / las_params["lpba"] + 0.1
    las_params["rlambda"] = 1.0 * las_params["mfbr"] / las_params["lpbr"] + 0.1

    to_print = ["I", "alim", "lpba", "mfba", "rlim", "lpbr", "mfbr", "ncurves0", "ncurves1"]
    sys.stderr.write("Using parameters %s\n" % " ".join(["%s:%s" % (key, las_params[key]) for key in to_print]))

    # Update parameters for makefb (which may depend on las parameters)
    update_existing(makefb_params, params)
    # also update "maxbits" which depends on "I"
    makefb_params["maxbits"] = params["I"] - 1
    makefb_params["out"] = las_params["fb"]

    makefb_cmd_line = [makefb]
    for (key, value) in makefb_params.items():
        makefb_cmd_line += ["-%s" % key, str(value)]
    # sys.stderr.write("Running: %s\n" % " ".join(makefb_cmd_line))
    subprocess.check_call(makefb_cmd_line)

    stats = LasStats()
    q0 = las_params["alim"]
    q_range = 1000
    q_inc = 0
    rels_wanted = int(primepi(las_params["lpba"]) + primepi(las_params["lpbr"]))
    # read the best time so far (if any)
    best_time = get_best_time ("las.best")
    if best_time >= 1e308:
       sys.stderr.write("Estimate %u relations needed\n" % rels_wanted)
    else:
       sys.stderr.write("Estimate %u relations needed, best time so far is %.0f seconds\n" % (rels_wanted, best_time))

    while stats.get_rels() < rels_wanted:
        # Set q0, q1
        outputfile = "las.%d.out" % q0
        las_params.update({"q0": q0, "q1": q0 + q_range, "out": outputfile})

        las_cmd_line = [las, "-allow-largesq"]
        for (key, value) in las_params.items():
            las_cmd_line += ["-%s" % key, str(value)]
        # sys.stderr.write("Running: %s\n" % " ".join(las_cmd_line))
        las_process = subprocess.Popen(las_cmd_line)
        (stdout, stderr) = las_process.communicate()
        if las_process.returncode != 0:
            raise Exception("las exited with status %d" % las_process.returncode)
        if not stats.parse_one_file(outputfile):
            raise Exception("Could not read statistics from %s" % outputfile)
        os.unlink(outputfile)
        cur_rels = stats.get_rels()
        if cur_rels == 0:
           cur_time = 0
        else:
           cur_time = stats.get_time(cur_rels)
        # sys.stderr.write("   Up to q=%u, estimate %u/%u relations in %.0f seconds\n" % (q0, cur_rels, rels_wanted, cur_time))

        # if the total time so far exceeds 1.5*best_time, we can stop here
        if cur_time > 1.5 * best_time:
           break

        # set q_inc so that we do about 10 sieving tests of length q_range
        if q_inc == 0:
           v = stats.relations_int.lastvalue[0]
           q_inc = int(rels_wanted / (10.0 * v))
        q0 += q_inc
    qmax = stats.get_qmax(rels_wanted)
    sievetime = stats.get_time(rels_wanted)
    if sievetime < best_time: # store the new best time
        set_best_time ("las.best", sievetime)
    sys.stderr.write("Estimated special-q up to %.0f and %.0f seconds for %d relations\n"
                      % (qmax, sievetime, rels_wanted))
    # Estimate how far the sieving should have gone to get the desired number
    # of relations
    return {'SIEVETIME': sievetime, 'RELATIONS': rels_wanted}

if __name__ == '__main__':
    param_file  = sys.argv[1]
    problem     = sys.argv[2]
    output_file = sys.argv[3]

    # Solve, gather measures and write to file.
    measures = run(param_file, problem)
    write_measures_to_file(output_file, measures)
