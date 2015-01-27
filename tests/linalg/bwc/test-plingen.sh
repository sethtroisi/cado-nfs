#!/usr/bin/env bash

set -e
set -x
# Create a fake sequence

# Note that if we arrive here, we are 64-bit only, since the GFP backends
# for bwc are explicitly disabled on i386 (for now -- most probably
# forever too).

wordsize=64

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

bindir="$1"
shift

dotest() {
    REFERENCE_SHA1="$1"
    shift

    : ${TMPDIR:=/tmp}
    TMPDIR=`mktemp -d $TMPDIR/plingen-test.XXXXXXXXXX`

    m="$1"; shift
    n="$1"; shift
    kmax="$1"; shift
    p="$1"; shift
    seed="$1"; shift

    unset thr
    unset ascii

    args=()
    mt_args=()
    for x in "$@" ; do
        case "$x" in
            plingen_program=*) eval "$x";;
            lingen_mpi_threshold*) mt_args=("${mt_args[@]}" "$x");;
            thr*) mt_args=("${mt_args[@]}" "$x"); thr="${x#thr=}";;
            *) args=("${args[@]}" "$x");
                if [[ "$x" =~ ascii ]] ; then
                    if [ "$p" -gt 1048576 ] ; then
                        echo "This test code support ascii tests only for small p" >&2
                        exit 1
                    fi
                    ascii=1
                fi
                ;;
        esac
    done

    # Hmmm, here's yet another dependency on bc -l...
    nbits_prime=$(echo "l($p)/l(2)+1" | bc -l | cut -d. -f1)
    nwords=$(echo "1+l($p)/l(2)/$wordsize"| bc -l)
    nwords=${nwords/.*/}

    : ${plingen_program:=plingen_p_$nwords}

    F="$TMPDIR/base"
    if [ "$ascii" ] ; then
        # The perl code below generates ascii test cases which are good
        # provided that p is small. Otherwise, the smallish coefficients
        # we generate are inappropriate and lead to failure, since
        # plingen guesses the length of the ascii input file.
        read -s -r -d '' code <<-'EOF'
            my ($m, $n, $kmax, $p, $seed) = @ARGV;
            my $u = int($seed / 1000);
            my $v = $seed % 1000;
            sub newx {
                    $u = ($u * 2) % 1048573;
                    $v = ($v * 3) % 1048573;
                    return ($u + $v) % $p;
                }

            for my $kmax (1..$kmax) {
                for my $i (1..$m) {
                    print join(" ", map { newx; } (1..$n)), "\n";
                }
                print "\n";
            }
            
EOF
        perl -e "$code" $m $n $((kmax/3)) $p $seed > $F
    else
        # generate $F with exactly ($kmax/3)*$m*$n*$nwords_per_gfp_elt
        # machine words of random data.
        read -s -r -d '' code <<-'EOF'
            my ($m, $n, $kmax, $nwords, $seed) = @ARGV;
            my $u = int($seed / 1000);
            my $v = $seed % 1000;
            sub newx {
                    $u = ($u * 2) % 1048573;
                    $v = ($v * 3) % 1048573;
                    return $u + $v;
                }
            for my $kmax (1..$kmax) {
                for my $i (1..$m) {
                    for my $j (1..$n) {
                        for my $s (1..$nwords) {
                            print pack("Q", newx);
                        }
                    }
                }
            }
            
EOF
        perl -e "$code" $m $n $((kmax/3)) $nwords $seed > $F
    fi

    G="$TMPDIR/seq.txt"
    cat $F $F $F > $G
    rm -f $F

    $bindir/linalg/bwc/$plingen_program m=$m n=$n prime=$p --afile $G "${args[@]}"
    [ -f "$G.gen" ]
    SHA1=$($SHA1BIN < $G.gen)
    SHA1="${SHA1%% *}"

    if [ "$REFERENCE_SHA1" ] ; then
        if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
            echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}" >&2
            exit 1
        fi
    else
        echo "========= $SHA1 ========"
    fi
    rm -f $G.gen

    mpi_bindir=$(perl -ne '/HAVE_MPI\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)

    if [ "$mpi_bindir" ] && [ "${mt_args[*]}" ] && [ "$thr" ] ; then
        set `echo $thr | tr 'x' ' '`
        if ! [ "$1" ] || ! [ "$2" ] ; then
            echo "Bad test configuration, thr should be of the form \d+x\d+ for MPI test" >&2
            exit 1
        fi
        njobs=$(($1*$2))
        $mpi_bindir/mpiexec -n $njobs $bindir/linalg/bwc/$plingen_program m=$m n=$n prime=$p --afile $G "${args[@]}" "${mt_args[@]}"

        [ -f "$G.gen" ]
        SHA1=$($SHA1BIN < $G.gen)
        SHA1="${SHA1%% *}"

        if [ "$REFERENCE_SHA1" ] ; then
            if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
                echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}" >&2
                exit 1
            fi
            rm -f $G.gen
        else
            echo "========= $SHA1 ========"
        fi
    fi

    rm -rf "$TMPDIR"
}

dotest "$@"

