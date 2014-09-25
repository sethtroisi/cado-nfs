from opal.core.io import *
import sys
import subprocess
import re
from math import log
from os import unlink

from report import LasStats

def update_existing(a, b):
    """ Update those keys in "a" which also exist in "b" with the values
    from "b"
    """
    for key in set(a) & set(b):
        a[key] = b[key]

def primepi(x):
    return x / (log(x) - 1)

def run(param_file, problem):
    "Run las with given parameters until the required number of relations is found."

    makefb = "/users/caramel/kruppaal/build/cado-nfs/normal/sieve/makefb"
    las = "/users/caramel/kruppaal/build/cado-nfs/normal/sieve/las"

    las_params = {
        "I": 11,
        "poly": "/tmp/work/test_run.polyselect2.poly",
        "fb": "/tmp/test_run.factorbase.roots.gz",
        "rlim": 50000,
        "alim": 100000,
        "lpbr": 22,
        "lpba": 22,
        "mfbr": 22,
        "mfba": 22,
        "rlambda": 1.2,
        "alambda": 1.2,
        "mt": 2
    }
    makefb_params = {
        "poly" : las_params["poly"],
        "alim": las_params["alim"],
        "maxbits": 10
    }

    params = read_params_from_file(param_file)

    # Update parameters for las
    update_existing(las_params, params)
    las_params["mfbr"] = max(las_params["mfbr"], las_params["lpbr"])
    las_params["mfba"] = max(las_params["mfba"], las_params["lpba"])
    
    to_print = ["I", "alim", "lpba", "mfba", "alambda", "rlim", "lpbr", "mfbr", "rlambda"]
    sys.stderr.write("Using parameters %s\n" % " ".join(["%s:%s" % (key, las_params[key]) for key in to_print]))

    # Update parameters for makefb (which may depend on las parameters)
    update_existing(makefb_params, params)
    makefb_params["out"] = las_params["fb"]

    makefb_cmd_line = [makefb]
    for (key, value) in makefb_params.items():
        makefb_cmd_line += ["-%s" % key, str(value)]
    # sys.stderr.write("Running: %s\n" % " ".join(makefb_cmd_line))
    subprocess.check_call(makefb_cmd_line)

    stats = LasStats()
    q0 = las_params["alim"]
    q_range = 1000
    q_inc = 10000
    rels_wanted = int(primepi(2**las_params["lpba"]) + primepi(2**las_params["lpbr"]))
    # sys.stderr.write("Estimate %f relations needed\n" % rels_wanted)
    
    while stats.get_rels() < rels_wanted:
        # Set q0, q1
        outputfile = "las.%d.out" % q0
        las_params.update({"q0": q0, "q1": q0 + q_range, "out": outputfile})

        las_cmd_line = [las]
        for (key, value) in las_params.items():
            las_cmd_line += ["-%s" % key, str(value)]
        # sys.stderr.write("Running: %s\n" % " ".join(las_cmd_line))
        las_process = subprocess.Popen(las_cmd_line)
        (stdout, stderr) = las_process.communicate()
        if las_process.returncode != 0:
            raise Exception("las exited with status %d" % las_process.returncode)
        if not stats.parse_one_file(outputfile):
            raise Exception("Could not read statistics from %s" % outputfile)
        unlink(outputfile)
        q0 += q_inc
    qmax = stats.get_qmax(rels_wanted)
    sievetime = stats.get_time(rels_wanted)
    sys.stderr.write("Estimated sq up to %f and %f seconds for %d relations\n"
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
