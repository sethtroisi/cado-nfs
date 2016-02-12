#!/bin/bash
if [ -f local.sh ]; then
	length=$( wc -l local.sh | cut -d " " -f 1 )
	if [ $length -ne 2 ]; then
		echo "Save your local.sh in local.sh.old."
		sleep 0.5
		cp local.sh local.sh.old
	fi
fi
echo 'CFLAGS="-O3 -DNDEBUG"' > local.sh
echo 'CXXFLAGS="-O3 -DNDEBUG"' >> local.sh
echo "Compile nfs-hd part of cado-nfs."
sleep 0.5
make makefb-hd special-q_sieve
cp parameters/nfs-hd/p6bd40.poly build/`hostname`/nfs-hd/.
cd build/`hostname`/nfs-hd
echo
echo "Launch makefb-hd (less than 10 minutes on a single 2GHz processor) for the p6bd40."
sleep 0.5
./makefb-hd -t 3 -fbb 1482911,1482911 -poly p6bd40.poly -lpb 25,25 -out sieve_base_p6bd40
echo "Launch special-q_sieve (less than 3 days on a single 2GHz processor) for the p6bd40."
sleep 0.5
./special-q_sieve -H 7,7,7 -fbb 1482911,1482911 -thresh 55,55 -poly p6bd40.poly -lpb 25,25 -q_side 1 -fb sieve_base_p6bd40 -q_range 1482937,2900828 -out sieve_1482937_2900828 -err /dev/null -gal 6
echo
echo "Check the relations"
sleep 0.5
cd ../../../../cado-nfs/nfs-hd/check_relations/
make build_tree=../../build/`hostname`/nfs-hd
cd ../../build/`hostname`/nfs-hd/
grep -v "#" sieve_1482937_2900828 > relations.raw
sort relations.raw | uniq > relations.uniq
./check_relations relations.uniq relations.true p6bd40.poly 25,25 relations.err
