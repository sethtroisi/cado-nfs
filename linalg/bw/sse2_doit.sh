#!/bin/zsh

D=~/Local/testmat

# MODULUS=9903520314283042199192993767
# A good FFT modulus (859*2^118+1, 128 bits).
MODULUS=17
# MODULUS=285451712094810683706092566195203997697
FFT_THRESHOLD=32
MULTISOLS=1

# 40/8/8/8 -> big pi without solution found.
MSIZE=2000
SRC=~/NDL/linalg
eval `make -s -C $SRC variables`
B="$BINARY_DIR"
if [ "$B" = "" ] ; then
	B=`pwd`/
fi
M=8
N=8
DENS=0.01

action() {
	echo "$@"
	"$@"
}

if [ -f "$1" ] ; then
	FILE=$1
	SPECIAL=$(head -c 2048 $FILE | md5sum | cut -c1-8)
	D="${D}.$SPECIAL"
	# REMOVE_D=yes
fi

rm -rf $D
mkdir $D

if [ "$1" = tiny ] ; then
	MSIZE=10
	DENS=0.6
	M=8
	N=8
fi

if [ -f "$1" ] ; then
	# works because we're zsh.
	grep ROWS "$1" | head -1 | tr -d ';' | read x x x MSIZE x x MODULUS
	cp "$1" $D/matrix.txt
	# unset FFT_THRESHOLD
else
	echo ${B}bw-random $MSIZE $MODULUS $DENS
	${B}bw-random $MSIZE $MODULUS $DENS > $D/matrix.txt
fi

M1=`expr $M - 1`
N1=`expr $N - 1`

action ${B}bw-balance --subdir $D
action ${B}bw-secure --subdir $D
action ${B}bw-prep --subdir $D $M $N

if [ "$1" = tiny ] ; then
	action ${B}bw-printmagma --subdir $D > $D/m.m
fi


for mi in {0..`expr $N1 / 8`} ; do
	i=`expr $mi \* 8`
	ni1=`expr $i + 7`
	if [ ! -f "$D/A-00-$(printf '%02d' $i)" ] ; then
		action ${B}bw-slave-mt --nthreads 2 --task slave --subdir $D 0,$i $M1,$ni1
	else
		echo "Column $i already done"
	fi
done

if [ "$FFT_THRESHOLD" != "" ] ; then
	action ${B}bw-master2 --subdir $D -t $FFT_THRESHOLD matrix.txt $M $N  | tee "$D/master.log"
else
	action ${B}bw-master-old --subdir $D matrix.txt $M $N  | tee "$D/master.log"
fi

# ./wrap_old_master_code.sh ${B}bw-master-old $D $MODULUS $MSIZE $M $N  | tee "$D/master.log"

# action ${B}bw-master --subdir $D | tee "$D/master.log"

ALLSOLS=`grep 'LOOK' "$D/master.log" | tr -d '[]' | cut -d\  -f6-`

#if [ -n $MULTISOLS ] ; then
	for X in `echo $ALLSOLS` ; do
		for mi in {0..`expr $N1 / 8`} ; do
			i=`expr $mi \* 8`
			ni1=`expr $i + 7`
			action ${B}bw-slave-mt --nthreads 2 --task mksol --subdir $D --sc $X 0,$i $M1,$ni1
		done
		action ${B}bw-gather --subdir $D $X --nbys 8
		cp "$D/W0$X" `dirname $FILE`
		cp "$D/MW0$X" `dirname $FILE`
	done
#else
#X=$(grep 'LOOK' "$D/master.log" | tail -1 | awk '// { print $6; }')
#
#for mi in {0..`expr $N1 / 8`} ; do
#	i=`expr $mi \* 8`
#	ni1=`expr $i + 7`
#        action ${B}bw-slave-mt --nthreads 2 --task mksol --subdir $D --sc $X 0,$i $M1,$ni1
#done
#
#action ${B}bw-gather --subdir $D $X --nbys 8
#fi

if [ "$REMOVE_D" = yes ] ; then
	F=`echo $FILE | sed -e s/matrix/solution/g`
	if [ "$F" = "$FILE" ] ; then
		F="${FILE}.solution"
	fi
	cp "$D/W"* `dirname $FILE`
	cp "$D/MW"* `dirname $FILE`
	cp "$D/W0$X" $F && rm -rf "$D"
fi
