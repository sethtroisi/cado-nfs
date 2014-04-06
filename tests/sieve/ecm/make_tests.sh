function make_factor_test {
  METHOD="$1"
  SIGMA="$2"
  SIGMANAME="_$SIGMA"
  if [ "$METHOD" = 1 ]
  then
    SIGMANAME=""
  fi
  B1="$3"
  B2="$4"
  declare -a METHODNAMES=(pm1 pp1 ecm ecmm12 ecmm16)
  METHODNAME="${METHODNAMES[$METHOD]}"
  if [ "$METHOD" = 1 ]
  then
    SIGMANAME=""
    METHODNAME="${METHODNAME}_${SIGMA/\//}"
  fi
  echo "find_test_combinations($B1,$B2,$SIGMA,$METHOD)" | gp -q ../../../sieve/ecm/ecm_avg_exp.gp >| "test_factor_${METHODNAME}${SIGMANAME}_${B1}_${B2}.inp"
}

function make_order_test {
  METHOD="$1"
  SIGMA="$2"
  SIGMANAME="_$SIGMA"
  if [ "$METHOD" = 1 ]
  then
    SIGMANAME=""
  fi
  PMIN="$3"
  PMAX="$4"
  declare -a METHODNAMES=(pm1 pp1 ecm ecmm12 ecmm16)
  METHODNAME="${METHODNAMES[$METHOD]}"
  if [ "$METHOD" = 1 ]
  then
    SIGMANAME=""
    METHODNAME="${METHODNAME}_${SIGMA/\//}"
  fi
  echo "print_orders($PMIN,$PMAX,$SIGMA,$METHOD)" | gp -q ../../../sieve/ecm/ecm_avg_exp.gp >| "test_order_${METHODNAME}${SIGMANAME}_${PMIN}_${PMAX}.inp"
}

make_factor_test 0 2 100 1000
make_factor_test 1 2/7 100 1000
make_factor_test 1 6/5 100 1000
make_factor_test 2 10 100 1000
make_factor_test 2 11 100 1000
make_factor_test 3 2 100 1000
make_factor_test 3 4 100 1000
make_factor_test 4 1 100 1000

for I in 10 11 12; do
  make_order_test 2 $I 1000 2000
  make_order_test 2 $I 1000000 1001000
  make_order_test 2 $I 1000000000 1000001000
done

for I in 4 5 6; do
  make_order_test 3 $I 1000 2000
  make_order_test 3 $I 1000000 1001000
done

make_order_test 4 1 1000 2000
make_order_test 4 1 1000000 1001000
