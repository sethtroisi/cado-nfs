#!/usr/bin/env bash
# print cpu time of cadofactor 
#
# Usage: $ cd "wdir"
#		 $ cpu_time [OPTIONS]
# OPTIONS: p: only polyselect
#		   s: only sieve (long time of computation)
#                  f: only filtering
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
  then  echo -n "CPU time for polyselect:    "
    if [ ! -f ${name}.polsel_out.* ] 2> /dev/null
      then echo "polynomial files were not found"
      else
         grep phase ${name}.polsel_out.* | sed "s/^.*phase took \([^s]*\).*$/+\1/g" | tr "\n" " " | cut -c2- | bc | f
         echo -n "   # rootsieve time:        "
         grep -h rootsieve ${name}.polsel_out.* | sed "s/^.*rootsieve took \([^s]*\).*$/+\1/g" | tr "\n" " " | cut -c2- | bc | f
         echo -n "   # of polyselect files:   "
         ls ${name}.polsel_out.* | wc -l
         echo -n "   # of found polynomials:  "
         grep tried ${name}.polsel_out.* | sed "s/^.*found \([^p]*\).*$/+\1/g" | tr "\n" " " | cut -c2- | bc
         echo -n "   # below maxnorm:         "
         grep tried ${name}.polsel_out.* | sed "s/^.*, \([0-9]*\) below.*$/+\1/g" | tr "\n" " " | cut -c2- | bc
    fi
fi

echo -n "   "
grep MurphyE ${name}.poly


# Sieve
if [[ -z $1 || $(expr $1 : '.*[s].*') != 0 ]]
  then echo -n "CPU time for sieve:         "
    if [ -f ${name}.rels ]
      then file=${name}.rels
      else if [ ! -f ${name}.rels.*.gz ] 2> /dev/null
        then echo "relations files were not found"
        else file=$(ls -f ${name}.rels.*.gz)
        fi
    fi
    if [ ! -z "$file" ]
      then zgrep "time" $file | sed "s/^.*time \([^s]*\).*$/+\1/g" | tr "\n" " " | cut -c2- | bc | f
    fi
    first=`ls ${name}.rels.*.gz | sed "s/^.*rels.\([0-9]*\)\-[0-9]*.gz$/\1/g" | sort -n | head -1`
    last=`ls ${name}.rels.*.gz | sed "s/^.*.rels.[0-9]*-\([0-9]*\).gz$/\1/g" | sort -n | tail -1`
    echo -n "   special-q range:         ${first}-${last}"
    lpba=`grep lpba ${name}.param | cut -c 6-`
    percent=`echo ${last}00/2^${lpba} | bc`
    echo " (${percent}% of LP bound)"
    echo -n "   # of sieve files:        "
    ls ${name}.rels.*.gz | wc -l
fi


echo -n "   # of relations found:    "
cat ${name}.nrels

# Filtering
if [[ -z $1 || $(expr $1 : '.*[f].*') != 0 ]]
  then echo -n "CPU time for dup1:      "
    if [ ! -f ${name}.dup1.log ] 2> /dev/null
      then echo "dup1 file was not found"
      else grep split ${name}.dup1.log | tail -1 | sed "s/^.*relations in \([^s]*\).*$/+\1/g" | tr "\n" " " | cut -c2- | bc | f
    fi
    echo -n "CPU time for dup2:      "
    if [ ! -f ${name}.dup2_*.log ] 2> /dev/null
      then echo "dup2 files were not found"
      else awk '/MB/ {x=$8;} /remaining/ {y+=x;} END {print y}' ${name}.dup2_*.log | f
    fi
    echo -n "   # of unique relations: "
    grep -h remaining ${name}.dup2*.log | sed "s/^At the end: \([0-9]*\) remaining.*/+\1/g" | tr "\n" " " | cut -c1- | sed "s/^+//g" | bc
    echo -n "CPU time for purge:     "
    if [ ! -f ${name}.purge.log ] 2> /dev/null
      then echo "purge file was not found"
      else grep MB ${name}.purge.log | tail -1 | sed "s/^.*relations in \([^s]*\).*$/+\1/g" | tr "\n" " " | cut -c2- | bc | f
    fi
    echo -n "CPU time for merge:     "
    if [ ! -f ${name}.merge.log ] 2> /dev/null
      then echo "merge file was not found"
      else grep "Total merge time" ${name}.merge.log | sed "s/^.*time: \([^s]*\).*$/\1/g" | f
    fi
    echo -n "CPU time for replay:    "
    if [ ! -f ${name}.replay.log ] 2> /dev/null
      then echo "replay file was not found"
      else grep "Total usage" ${name}.replay.log | sed "s/^.*time \([^s]*\).*$/\1/g" | f
    fi
fi

# Linear Algebra
if [[ -z $1 || $(expr $1 : '.*[l].*') != 0 ]]
  then echo -n "Real time for linear algebra: "
  if [ ! -f ${name}.bwc.log ] 2> /dev/null
    then echo "linear algebra logfile were not found"
    else line=$(grep -n Done. ${name}.bwc.log | head -1 | cut -d':' -f1)
    krylov=$(head -$line ${name}.bwc.log | grep 's/iter' | grep krylov | grep "000 " | sed "s/^.* \[\(.*\) s\/it.*$/+1000*\1/g" | tr "\n" " " | cut -c2- | bc | f)
    lingen=$(grep "^Total computation took" ${name}.bwc.log | cut -d' ' -f4 | f)
    mksol=$(grep 's/iter' ${name}.bwc.log | grep mksol | grep "000 " | sed "s/^.* \[\(.*\) s\/it.*$/+1000*\1/g" | tr "\n" " " | cut -c2- | bc | f)
    echo "$krylov - $lingen - $mksol (krylov - lingen - mksol)"
  fi	
fi


# Sqrt
if [[ -z $1 || $(expr $1 : '.*[r].*') != 0 ]]
  then echo -n "CPU time for sqrt (one root): "
    if [ ! -f ${name}.sqrt.log ] 2> /dev/null
      then echo "sqrt logfile were not found"
    else
      rat=$(grep Rational ${name}.sqrt.log | head -1 | sed "s/.*: \(.*\)ms$/\1/g")
      if [ ! -z $rat ]
      then rat2=$(echo "scale=3; $rat/1000" | bc -l | f)
           alg=$(grep Algebraic ${name}.sqrt.log | head -1 | cut -d' ' -f6 | f)
           echo "${rat2} - ${alg}  (ratsqrt - algsqrt)"
      else echo "sqrt not calculated"
      fi
    fi
fi
