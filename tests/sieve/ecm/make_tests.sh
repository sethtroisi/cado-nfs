function make_test {
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
  echo "find_test_combinations($B1,$B2,$SIGMA,$METHOD)" | gp -q ../../../sieve/ecm/ecm_avg_exp.gp >| "${METHODNAME}${SIGMANAME}_${B1}_${B2}_test.inp"
}

make_test 0 2 100 1000
make_test 1 2/7 100 1000
make_test 1 6/5 100 1000
make_test 2 10 100 1000
make_test 2 11 100 1000
make_test 3 2 100 1000
make_test 3 4 100 1000
make_test 4 1 100 1000
