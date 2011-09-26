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

name=$(grep '^name=' *.param | cut -d"=" -f2)
#name=$2

function f
{
perl -e 'sub format_dhms {
  my $sec = shift;
  my ($d, $h, $m);
  $d = int ( $sec / 86400 ); $sec = $sec % 86400;
  $h = int ($sec / 3600 ); $sec = $sec % 3600;
  $m = int ($sec / 60 ); $sec = $sec % 60;
  if ($d!=0) {
    return "${d}d${h}h${m}m${sec}s";
  } else {
    if ($h!=0) {
      return "${h}h${m}m${sec}s";
    } else {
      if ($m!=0) {
        return "${m}m${sec}s";
      } else {
        return "${sec}s";
      }
    }
  }
}
my $var=<STDIN>; print format_dhms($var)."\n"'
}


# Polyselect
if [[ -z $1 || $(expr $1 : '.*[p].*') != 0 ]]
  then  echo -n "CPU time for polyselect:      "
    if [ ! -f ${name}.kjout.* ] 2> /dev/null
      then echo "polynomial files were not found"
      else grep phase ${name}.kjout.* | sed "s/^.*phase took \(\S*\)s.*$/\1/g" | tr "\n" "+" | sed "s/^\(.*\)+$/\1\n/" | bc | f
    fi
fi


# Sieve
if [[ -z $1 || $(expr $1 : '.*[s].*') != 0 ]]
  then echo -n "CPU time for sieve:           "
    if [ -f ${name}.rels ]
      then file=${name}.rels
      else if [ ! -f ${name}.rels.*.gz ] 2> /dev/null
        then echo "relations files were not found"
        else file=$(ls -f ${name}.rels.*.gz)
        fi
    fi
    if [ ! -z "$file" ]
      then zgrep "time" $file | sed "s/^.*time \(\S*\)s.*$/\1/g" | tr "\n" "+" | sed "s/^\(.*\)+$/\1\n/" | bc | f
    fi
fi


# Linear Algebra
if [[ -z $1 || $(expr $1 : '.*[l].*') != 0 ]]
  then echo -n "Real time for linear algebra: "
  if [ ! -f ${name}.bwc.log ] 2> /dev/null
    then echo "linear algebra logfile were not found"
    else line=$(grep -n Done. ${name}.bwc.log | head -1 | cut -d':' -f1)
    krylov=$(head -$line ${name}.bwc.log | grep 's/iter' | grep krylov | grep "000 " | sed "s/^.* \[\(.*\) s\/iter\]$/\1/g" | tr "\n" "+" | sed "s/^\(.*\)+$/(\1)*1000\n/" | bc | f)
    lingen=$(grep "^Total computation took" ${name}.bwc.log | cut -d' ' -f4 | f)
    mksol=$(grep 's/iter' ${name}.bwc.log | grep mksol | grep "000 " | sed "s/^.* \[\(.*\) s\/iter\]$/\1/g" | tr "\n" "+" | sed "s/^\(.*\)+$/(\1)*1000\n/" | bc | f)
    echo "$krylov - $lingen - $mksol (krylov - lingen - mksol)"
  fi	
fi


# Sqrt
if [[ -z $1 || $(expr $1 : '.*[r].*') != 0 ]]
  then echo -n "CPU time for sqrt (one root): "
    if [ ! -f ${name}.sqrt.log ] 2> /dev/null
      then echo "sqrt logfile were not found"
    else
      rat=$(grep Rational ${name}.sqrt.log | head -1 | sed "s/.* at \(\S*\)ms$/\1/g")
      if [ ! -z $rat ]
      then rat2=$(echo "scale=3; $rat/1000" | bc -l | f)
           alg=$(grep Algebraic ${name}.sqrt.log | head -1 | cut -d' ' -f6 | f)
           echo "${rat2} - ${alg}  (ratsqrt - algsqrt)"
      else echo "sqrt not calculated"
      fi
    fi
fi
