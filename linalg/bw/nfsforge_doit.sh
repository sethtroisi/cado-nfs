#!/bin/zsh

D=~/Local/testmat

# MODULUS=9903520314283042199192993767
# A good FFT modulus (859*2^118+1, 128 bits).
# MODULUS=285451712094810683706092566195203997697
# FFT_THRESHOLD=512

SRC=~/NDL/linalg
eval `make -s -C $SRC variables`
B="$BINARY_DIR"
if [ "$B" = "" ] ; then
	B=`pwd`/
fi
M=2
N=2
DENS=0.002

action() {
	echo "$@"
	"$@"
}

if ! [ -f "$1" ] ; then
	echo "need a file"
	exit 1
fi
FILE="$1"
SPECIAL=$(head -c 2048 $FILE | md5sum | cut -c1-8)
D="${D}.$SPECIAL"
# REMOVE_D=yes

rm -rf $D
mkdir $D

M=16
N=16
# works because we're zsh.
grep ROWS "$1" | head -1 | tr -d ';' | read x x x MSIZE x x MODULUS
cp "$1" $D/matrix.txt
unset FFT_THRESHOLD

M1=`expr $M - 1`
N1=`expr $N - 1`

action ${B}bw-balance --subdir $D
action ${B}bw-secure --subdir $D
action ${B}bw-prep --subdir $D $M $N

if [ "$1" = tiny ] ; then
	${B}bw-printmagma --subdir $D > $D/m.m
fi


for i in {0..$N1} ; do
	if [ ! -f "$D/A-00-$(printf '%02d' $i)" ] ; then
		action ${B}bw-slave-mt --nthreads 2 --task slave --subdir $D $i
	else
		echo "Column $i already done"
	fi
done

if [ "$FFT_THRESHOLD" != "" ] ; then
	action ${B}../../matlingen/bw-master2 --subdir $D -t $FFT_THRESHOLD matrix.txt $M $N  | tee "$D/master.log"
else
	action ${B}../../matlingen/bw-master --subdir $D matrix.txt $M $N  | tee "$D/master.log"
fi

# ./wrap_old_master_code.sh ${B}bw-master-old $D $MODULUS $MSIZE $M $N  | tee "$D/master.log"

# action ${B}bw-master --subdir $D | tee "$D/master.log"

X=$(grep 'LOOK' "$D/master.log" | tail -1 | awk '// { print $6; }')

for i in {0..$N1} ; do
        action ${B}bw-slave-mt --nthreads 2 --task mksol --subdir $D --sc $X $i
done

${B}bw-gather --subdir $D $X

if [ "$REMOVE_D" = yes ] ; then
	F=`echo $FILE | sed -e s/matrix/solution/g`
	if [ "$F" = "$FILE" ] ; then
		F="${FILE}.solution"
	fi
	cp "$D/W0$X" $F && rm -rf "$D"
fi
