#!/bin/zsh

if [ -x /usr/local/sbin/tell_cpu ] ; then
	/usr/local/sbin/tell_cpu -q max
fi

D=~/Local/testmat2

M=9903520314283042199192993767
S=1000
make variables
eval `make variables`
B="$BINARY_DIR"
MN=2
DENS=0.05

action() {
	echo "$@"
	"$@"
}

new_master() {
	echo -n "new master code:\t"
	time ${B}bw-master --subdir $D > "$D/master.log"
}

old_master() {
	echo -n "old master code:\t"
	(cd "$D" ; time ${B}bw-master-old 0 > "$D/master-old.log")
}

old_master
new_master


X=$(grep 'LOOK' "$D/master.log" | tail -1 | awk '// { print $6; }')

J1=`expr $MN - 1`

finish() {
for i in {0..$J1} ; do
        action ${B}bw-slave-mt --nthreads 2 --task mksol --subdir $D --sc $X $i
done

${B}bw-gather --subdir $D $X
}

echo -n "Checking the solution: "

finish 2>&1  | grep -q 'is a solution' && echo ok

if [ -x /usr/local/sbin/tell_cpu ] ; then
	/usr/local/sbin/tell_cpu -q min
	/usr/local/sbin/tell_cpu -q dyn
fi

