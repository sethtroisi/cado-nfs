#!/bin/zsh

BIN=$1
D=$2
MODULUS=$3
MSIZE=$4

M=$5
N=$5

action() {
	echo "$@"
	"$@"
}

M1=`expr $M - 1`
N1=`expr $N - 1`
MN1=`expr $M + $N - 1`

gcc -Wall -O -o /tmp/dec2bin dec2bin.c -lgmp
gcc -Wall -O -o /tmp/bin2dec dec2bin.c -lgmp

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

action "$BIN" 0

for j in {00..0$MN1} ; do
	/tmp/bin2dec 128 rF-{00..$M1}-$j.* | tac > F$j
done
