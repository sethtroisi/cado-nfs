#!/usr/bin/env python3
#
# Example for the p59 in tests/test_full_p59
# scripts/descent_finish.py
# ...

import os
import subprocess
import sys
import argparse
import re
import math
import tempfile
import shutil

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Descent final step for DLP")
    # Required
    parser.add_argument("--lpb", help="Large prime bound", required=True,
            type=int)
    parser.add_argument("--cadobindir", help="Cado build directory",
            required=True, type=str)
    parser.add_argument("--rels", help="File of relations of the descent step",
            required=True, type=str)
    parser.add_argument("--renumber", help="Renumber file",
            required=True, type=str)
    parser.add_argument("--debug_ren", help="File given by debug_renumber",
            required=True, type=str)
    parser.add_argument("--poly", help="Polynomial file",
            required=True, type=str)
    parser.add_argument("--badidealinfo", help="Badideal info file",
            required=True, type=str)
    parser.add_argument("--gorder", help="Group order (a.k.a. ell)",
            required=True, type=int)
    parser.add_argument("--nsm0", help="Number of SM on side 0", required=True,
            type=int)
    parser.add_argument("--nsm1", help="Number of SM on side 1", required=True,
            type=int)
    parser.add_argument("--smexp0", help="SM exponant on side 0", required=True,
            type=int)
    parser.add_argument("--smexp1", help="SM exponant on side 1", required=True,
            type=int)
    parser.add_argument("--log", help="File with known logs", required=True,
            type=str)
    # Optional
    parser.add_argument("--workdir", help="Temporary working directory")
    parser.add_argument("--mt", help="Multi-thread level of reconstructlog",
            default=4)

    args = parser.parse_args()
    wdir = args.workdir
    if wdir == None:
        wdir = tempfile.mkdtemp(dir="/tmp")
    print ("Working directory is: " + wdir)
    descrels_filename = args.rels
    LPB = 1<<args.lpb

    # Read descent relations
    descrels = []
    with open(descrels_filename, 'r') as relfile:
        for line in relfile:
            if line[0] == '#':
                continue;
            else:
                descrels.append(line.rstrip())
    nrels = len(descrels)
    print ("Have read " + str(nrels) + " relations.")

    # Create fake relations: remove large primes
    fakerels = []
    extraprimes = []
    for rel in descrels:
        r = rel.split(':')
        ss = [ None, None]
        extra = [ None, None]
        for side in range(2):
            rr = r[side+1].split(',')
            ss[side] = []
            extra[side] = []
            for p in rr:
                if int(p, 16) < LPB:
                    ss[side].append(p)
                else:
                    extra[side].append(p)
        fr = r[0]
        for side in range(2):
            fr = fr + ":"
            for i in range(len(ss[side])):
                fr = fr + ss[side][i]
                if i < len(ss[side])-1:
                    fr = fr + ","
        fakerels.append(fr)
        extraprimes.append(extra)
    setextraprimes = set()
    for extra in extraprimes:
        for side in range(2):
            for p in extra[side]:
                setextraprimes.add(p)
    listextraprimes = list(setextraprimes)
    listextraprimes.sort()
    print ("They include " + str(len(listextraprimes)) + " primes above lpb")
    fakerels_filename = wdir + "/fakerels"
    with open(fakerels_filename, 'w') as file:
        for rel in fakerels:
            file.write(rel + '\n')

    # Call dup2 to renumber those relations to get them renumbered
    # It works in place, so we do a copy first.
    fakerels_ren_filename = wdir + "/fakerels_ren"
    subprocess.check_call(["cp", fakerels_filename,
        fakerels_ren_filename])
    subprocess.check_call([args.cadobindir + "/filter/dup2",
        "-dl",
        "-nrels", str(nrels),
        "-renumber", args.renumber,
        "-badidealinfo", args.badidealinfo,
        fakerels_ren_filename
        ])

    # In the renumber file, the line number is not exactly the index.
    # Let us rely on debug_renumber's output to find the last used index
    # in the renumber table.
    last_renumber_line = None
    with subprocess.Popen(["tail", "-10", args.debug_ren],
            stdout=subprocess.PIPE) as p:
        for line in p.stdout.readlines():
            ll = line.decode("utf-8")
            if ll[0] == 'i':
                last_renumber_line = ll.rstrip()
    last_index = int(last_renumber_line.split()[0].split('=')[1], 16)

    # Now, we can assign fresh indices to new primes
    # We fix the relations and the create a new renumber table
    # accordingly.
    dictextraprimes = {}
    for i in range(len(listextraprimes)):
        dictextraprimes[listextraprimes[i]] = last_index+1+i
    descrels_ren_filename = wdir + "/descrels_ren"
    # relations
    i = 0
    with open(fakerels_ren_filename, 'r') as infile:
        with open(descrels_ren_filename, 'w') as outfile:
            # add first line like in purged.gz, just to give the number
            # of relations to reconstructlog.
            outfile.write("# " + str(nrels) + "0 0\n")
            for rel in infile:
                outfile.write(rel.rstrip())
                for side in range(2):
                    for p in extraprimes[i][side]:
                        ind = hex(dictextraprimes[p])[2:]
                        outfile.write("," + ind)
                outfile.write("\n")
                i = i+1
    # renumber file
    largest_prime = int(listextraprimes[len(listextraprimes)-1], 16)
    newlpb = math.ceil(math.log(largest_prime, 2))
    newrenumber_filename = wdir + "/new_renumber"
    with open(newrenumber_filename, 'w') as file: 
        subprocess.check_call(["gzip", "-cd", args.renumber], stdout=file)
    with open(newrenumber_filename, 'r+') as file:
        line1 = file.readline()
        newline1 = line1;
        # The end of the line looks like " 26 26\n", where 26 is the lpb
        ll = len(line1)
        assert newline1[ll-1] == '\n'
        assert newline1[ll-4] == ' '
        assert newline1[ll-7] == ' '
        nlpb = str(newlpb)
        assert len(nlpb) == 2
        newline1 = newline1[0:ll-6] + nlpb + " " + nlpb + "\n"
        assert len(newline1) == ll
        file.seek(0)
        file.write(newline1)
    with open(newrenumber_filename, 'a') as file:
        for p in listextraprimes:
            pp = hex(int(p, 16)+1)[2:]
            file.write(pp+'\n')

    # Prepare the call to reconstructlog
    # Need to modifiy the last lines of known logs that contain the
    # virtual logs of the SMs
    nsm = args.nsm0 + args.nsm1
    newlog_filename = wdir + "/new_log"
    count_sm=0
    with open(args.log, 'r') as infile:
        with open(newlog_filename, 'w') as outfile:
            for line in infile:
                if not "SM" in line:
                    outfile.write(line)
                else:
                    count_sm = count_sm + 1
                    ll = line.split()
                    assert len(ll) == 5
                    ind = int(ll[0], 16)
                    ind = ind + len(listextraprimes)
                    ll[0] = hex(ind)[2:]
                    outfile.write(ll[0] + " " + ll[1] + " "
                            + ll[2] + " " + ll[3] + " " + ll[4] + "\n")
    assert nsm == count_sm

    result_filename = wdir + "/result_log"
    fake_relsdel_filename = wdir + "/fake_relsdel"
    subprocess.check_call(["touch", fake_relsdel_filename])
    subprocess.check_call([args.cadobindir + "/filter/reconstructlog-dl",
        "-gorder", str(args.gorder),
        "-sm0", str(args.nsm0),
        "-sm1", str(args.nsm1),
        "-smexp0", str(args.smexp0),
        "-smexp1", str(args.smexp1),
        "-nrels", str(nrels),
        "-renumber", newrenumber_filename,
        "-log", newlog_filename,
        "-logformat", "reconstruct",
        "-out", result_filename,
        "-poly", args.poly,
        "-purged", descrels_ren_filename,
        "-relsdel", fake_relsdel_filename,
        "-mt", str(args.mt)])

    print ("Workdir is " + wdir)
    print ("Result is in " + result_filename)
