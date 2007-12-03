#!/bin/zsh

D=~/Local/testmat2

MODULUS=9903520314283042199192993767
MSIZE=1000
make variables
eval `make variables`
B="$BINARY_DIR"
if [ "$B" = "" ] ; then
	B=`pwd`/
fi
M=2
N=2
DENS=0.05

action() {
	echo "$@"
	"$@"
}

rm -rf $D
mkdir $D
${B}bw-random $MSIZE $MODULUS $DENS > $D/matrix.txt
action ${B}bw-balance --subdir $D
action ${B}bw-secure --subdir $D
action ${B}bw-prep --subdir $D $M $N

M1=`expr $M - 1`
N1=`expr $N - 1`

for i in {00..0$N1} ; do
	action ${B}bw-slave-mt --nthreads 2 --task slave --subdir $D $i
done

gcc -O -o /tmp/dec2bin dec2bin.c -lgmp
gcc -O -o /tmp/bin2dec dec2bin.c -lgmp

for i in {00..0$N1} ; do
	for j in {00..0$N1} ; do
		/tmp/dec2bin 128 < "$D/A-$i-$j" > "$D/A-$i-$j.000"
	done
done

NI=`expr $MSIZE / $M + $MSIZE / $N + 2 \* $M / $N + 2 \* $N / $M + 10`
LOTS=`expr 10 \* $NI`
echo ":$M:$N:$MSIZE:$NI:$LOTS:$MODULUS" > "$D/WORKING"

cd "$D"
ln -s . run
ln -s . output
ln -s . 0

