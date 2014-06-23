#!/usr/bin/env bash

# This script downloads, builds, and runs cado-nfs on several machines.
# The machines from the gcc compile farm are taken as an example. Tests
# are done in parallel on the different nodes. Once all jobs are
# finished, a digest is sent by mail, with a tar.gz attachment containing
# all the results. The output can probably be improved (e.g. give
# warnings in a text attachments, not buried inside a tar.gz).  The
# requirement is that cmake is present on the target machines (as per the
# companion script, here gcc-script.sh).

# The companion script gcc-script.sh is of course required.

# In order to access the machines through ssh, it is mandatory to have
# proper support in  ~/.ssh/config. For the gcc compile farm, e.g.:
# GET http://gcc.gnu.org/wiki/CompileFarm | perl ~thome/.ssh/parse-cfarm.pl

# results with cado-nfs-1.1.tar.gz
gccs=(\
gcc10 # ok
gcc11 # ok
gcc12 # ok
gcc13 # ok
gcc14 # ok
gcc15 # ok
gcc16 # ok
gcc17 # ok
gcc20 # ok
# gcc33 Connection refused
# gcc34 Connection refused
# gcc35 Connection refused
# gcc36 Connection refused
# gcc37 Connection refused
gcc38 # ok
# gcc40 Connection timed out
# gcc41 Connection refused
# gcc42 Connection refused
# gcc43 Connection refused
gcc45 # ok
gcc46 # ok
gcc47 # ok
# gcc50 Connection timed out
gcc51 # ok
# gcc52 Connection refused
# gcc53 Connection timed out
gcc54 # ok
# gcc55 Connection refused
# gcc56 Connection timed out
# gcc57 Connection timed out
gcc60 # ok
gcc61 # ok
# gcc62 Connection timed out
gcc63 # ok
gcc64 # ok with egcc 4.2.4 (GMP compiled with egcc too)
gcc66 # ok
gcc70 # coredump in dispatch -> bug in pthreads?
# gcc100 Connection timed out
# gcc101 Connection timed out
gcc110
# gcc200 Connection timed out
# gcc201 Connection timed out
)

tmpdir=`mktemp -d /tmp/XXXXXXXX`

check_paul_on_one_node() {
    host=$1
    ssh -x -oStrictHostKeyChecking=no $host ls /home/zimmerma/bin 2>&1 | sed -e "s/^/$host:/"
}

run_on_one_node() {
    host=$1
    ssh -R1234:localhost:22 $host < "`dirname $0`/gcc-script.sh" > $tmpdir/$host.output 2>&1
}

for host in "${gccs[@]}" ; do
    run_on_one_node $host & :
done
wait
for host in "${gccs[@]}" ; do
    tail -2 $tmpdir/$host.output | sed -e s,^,$host:,
    echo
done > $tmpdir/digest.txt

# TODO: Include a log of warnings and/or compilation errors.

TODAY=`date +%Y%m%d`
TAR=$tmpdir/nightlytest-$TODAY.tar.gz

(cd $tmpdir ; tar cz --transform s,^,nightlytest-$TODAY/, -f - *.output digest.txt > $TAR)

mutt -s "nightly test $TODAY" -i $tmpdir/digest.txt -a $TAR -- `id -u -n` < /dev/null >/dev/null 2>&1

rm -rf $tmpdir
