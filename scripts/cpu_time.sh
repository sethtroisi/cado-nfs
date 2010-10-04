#!/bin/bash
# print cpu time of cadofactor 
#
# Usage: $ cd "wdir"
#		 $ cpu_time [OPTIONS]
# OPTIONS: p: only polyselect
#		   s: only sieve (long time of computation)
#		   l: only linalg
#		   r: only sqrt
#	       ps: only polyselect and sieve
#		   ...

name=$(grep '^name=' *.param | sed "s/^name=\(\S*\)$/\1/g")

function f
{
perl -e 'sub format_dhms {
  my $sec = shift;
  my ($d, $h, $m);
  $d = int ( $sec / 86400 ); $sec = $sec % 86400;
  $h = int ($sec / 3600 ); $sec = $sec % 3600;
  $m = int ($sec / 60 ); $sec = $sec % 60;
  return "$d"."d:$h"."h:$m"."m:$sec"."s";
}
my $var=<STDIN>; print format_dhms($var)."\n"'
}

if [[ -z $1 || $(expr $1 : '.*[p].*') != 0 ]]
  then echo -n "CPU time for polyselect:      "
    if [ ! -f ${name}.kjout.* ] 2> /dev/null
      then stty -echo; read c; stty echo
        ls ${name}_$c.kjout.* | while read line; do new_line=$(echo $line | sed "s/^${name}_$c\.\(.*\)/${name}\.\1/g"); cp $line $new_line; done
    fi
    grep phase ${name}.kjout.* | sed "s/^.*phase took \(\S*\)s.*$/\1/g" | tr "\n" "+" | sed "s/^\(.*\)+$/\1\n/" | bc | f
fi

if [[ -z $1 || $(expr $1 : '.*[s].*') != 0 ]]
  then echo -n "CPU time for sieve:           "
    if [ -f ${name}.rels ]
      then b=${name}.rels
      else b=$(ls -f ${name}.rels.*.gz)
    fi
  zgrep "time" $b | sed "s/^.*time \(\S*\)s.*$/\1/g" | tr "\n" "+" | sed "s/^\(.*\)+$/\1\n/" | bc | f
fi

if [[ -z $1 || $(expr $1 : '.*[l].*') != 0 ]]
  then echo "CPU or real time for linear algebra"
    num_lines=$(grep -n Done. ${name}.bwc.stderr | sed "s/^\([0-9]*\):Done./\1-1/g")
    if [[ $(echo $num_lines | wc -w) == 2 ]]
      then echo -n "        krylov [real time]:   "
        num_line_v=$(echo $num_lines | cut -d' ' -f 1 | bc)
        line_v=$(head -$num_line_v ${name}.bwc.stderr | tail -1)
        v0=$(echo $line_v | sed "s/^.*ETA (N=\(\S*\)): .*/\1/g")
        v1=$(echo $line_v | sed "s/^.*\[\(\S*\) s\/iter\]$/\1/g")
        t=$(echo $v0*$v1 | bc | f)
        echo "$t ($v1 s/iter)"
        echo -n "        lingen:               "
        grep "^Total computation took" ${name}.bwc.stderr | sed "s/^.*took \(\S*\)$/\1/g" | f
        echo -n "        mksol [real time]:    "
        num_line_s=$(echo $num_lines | cut -d' ' -f 2 | bc)
        line_s=$(head -$num_line_s ${name}.bwc.stderr | tail -1)
        s0=$(echo $line_s | sed "s/^.*ETA (N=\(\S*\)): .*/\1/g")
        s1=$(echo $line_s | sed "s/^.*\[\(\S*\) s\/iter\]$/\1/g")
        t=$(echo $s0*$s1 | bc | f)
        echo "$t ($s1 s/iter)"
      else echo "    mksol not ended or linalg has been rerun with checkpoints"
    fi	
fi

if [[ -z $1 || $(expr $1 : '.*[r].*') != 0 ]]
  then echo "CPU time for sqrt (one root): "
    echo -n "        ab pairs:             "
    grep dependency ${name}.sqrt.stderr | head -1 | cut -d' ' -f 6 | f
    echo -n "        ratsqrt:              "
    rat=$(grep Rational ${name}.sqrt.stderr | head -1 | sed "s/.* at \(\S*\)ms$/\1/g")
    echo "scale=3; $rat/1000" | bc -l | f
    echo -n "        algsqrt:              "
    grep Algebraic ${name}.sqrt.stderr | head -1 | cut -d' ' -f6 | f
fi
