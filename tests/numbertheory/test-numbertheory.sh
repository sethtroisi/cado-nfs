#!/bin/bash

if [ "$CADO_DEBUG" ] ; then
    set -x
fi

nfields=10
primelist="2 3 5 7 11 13 17 19 23"

while [ $# -gt 0 ] ; do
    if [ "$1" = "--binary" ] ; then
        binary=$2
        shift
    elif [ "$1" = "--magma" ] ; then
        magma=$2
        shift
    elif [ "$1" = "--fields" ] ; then
        nfields=$2
        shift
    elif [ "$1" = "--primelist" ] ; then
        primelist=$2
        shift
    else
        echo "unexpected arg: $1" >&2
        exit 1
    fi
    shift
done

if ! [ "$binary" ] || ! [ -x "$binary" ] ; then
    echo "binary=$binary: not found or not executable" >&2
    exit 1
fi
if ! [ "$magma" ] || ! [ -x "$magma" ] ; then
    echo "magma=$magma: not found or not executable" >&2
    exit 1
fi

input_file_for_binary() {
    d=$1
    shift
    echo $d
    echo "$@"
}

read -s -r -d '' perl_code <<-'EOF'
    while (<>) {
        /maximal order is/ && do { $x=1; next; };
        s/^\s*//;
        print if $x;
        do { print "\n"; last; } if $x && /^$/;
    }
    while (<>) {
        s/^\s*//;
        print if /^[\d\s\/]*$/;
        /multiplicity of (\d+)/ && do { print "$1\n"; };
    }

EOF

magma_code() {
cat <<-'EOF'
    f:=Polynomial(StringToIntegerSequence(COEFFS));
    n:=Degree(f);
    K<alpha>:=NumberField(f);
    rr:=Split(RESULT," ");
    r_order:=rr[1..n^2];
    rr:=rr[n^2+1..#rr];
    assert #rr mod (n^2+1) eq 0;
    nideals:=#rr div (n^2+1);
    basis:=Matrix(Rationals(),n,n,[StringToRational(x):x in r_order]);
    ideal_bases:=[<Matrix(Rationals(),n,n,[StringToRational(x):x in
    rr[(k-1)*(n^2+1)+1..k*(n^2+1)-1]]),StringToInteger(rr[k*(n^2+1)])> : k in [1..nideals]];
    O:=Order([LeadingCoefficient(f)*alpha]);
    OK:=MaximalOrder(K);
    try
        testing:="p-maximal order";
        myO:=Order([FieldOfFractions(O)!Eltseq(r):r in Rows(basis)]);
        assert Gcd(Index(pMaximalOrder(myO,p),myO), p) eq 1;
        testing:="factorization of p";
        myO:=Order([FieldOfFractions(O)!Eltseq(r):r in Rows(basis)]);
        myIs:=[<ideal<myO|[FieldOfFractions(O)!Eltseq(r):r in
        Rows(x[1])]>, x[2]>:x in ideal_bases];
        myI_OKs:=[<ideal<OK|I[1]>,I[2]>:I in myIs];
        assert Seqset(Decomposition(OK,p)) eq Seqset(myI_OKs);
        printf "ok p:=%o; d:=%o; COEFFS:=\"%o\";\n", p, n, COEFFS, RESULT;
    catch e
        printf "FAILED[%o] p:=%o; d:=%o; COEFFS:=\"%o\"; RESULT:=\"%o\";\n", testing, p, n, COEFFS, RESULT;
    end try;
    quit;
EOF
}


computed_basis() {
    $binary "$@" | perl -e "$perl_code" | xargs echo
}

primes=($(echo $primelist))

set -e

for n in `seq 1 $nfields` ; do
    d=$((RANDOM % 10))
    let d+=2
    coeffs=()
    x=()
    for i in `seq 0 $d` ; do
        if [ $i -lt $d ] ; then
            x=("${x[@]}" $((RANDOM-10000)))
        else
            while true ; do
                lc=$RANDOM
                if [ $lc != 0 ] ; then break; fi
            done
            x=("${x[@]}" $lc)
        fi
    done
    set -- "${x[@]}"

    for p in "${primes[@]}" ; do
        res=$(computed_basis <(input_file_for_binary $d "${x[@]}") $p)
        set -- "${x[@]}"
        displayed=$($magma -b <(echo "COEFFS:=\"$*\"; RESULT:=\"$res\"; p:=$p;") <(magma_code)  < /dev/null)
        echo $displayed
        if echo $displayed | grep -q FAILED ; then
            exit 1
        fi
    done
done
