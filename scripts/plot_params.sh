# Run inside the parameter directory

PARAMDIR=.
if [ -n "$1" ]
then
  PARAMDIR="$1"
fi

if ! test -d "$PARAMDIR" || ! test -f "$PARAMDIR/params.c90"
then
  echo Please specify the parameter directory
  exit 1
fi

KEYWORDS=( degree qrange I polyselect.P polyselect.admax polyselect.adrange lim0 lim1 lpb0 lpb1 )
declare -a OUTPUTFILES

for S in "${KEYWORDS[@]}"
do
  OUTPUTFILE="$S.png"
  OUTPUTFILES=( "${OUTPUTFILES[@]}" "$OUTPUTFILE" )
  echo "Creating $OUTPUTFILE"
  grep -H "^[^#]*$S[[:space:]]*=" "$PARAMDIR"/params.c[1-9][0-9] "$PARAMDIR"/params.c[12][0-9][0-9] | sed "s/^.*params.c\([0-9]*\):.*$S[[:space:]]*=[[:space:]]*\([0-9a-zA-Z.e-]*\).*/\1 \2/;" >| /tmp/$S.data
  echo "set terminal png; set log y 2; set output \"$OUTPUTFILE\"; plot \"/tmp/$S.data\" " | gnuplot
done
echo "Creating montage.png"
montage -geometry 640x480 "${OUTPUTFILES[@]}" montage.png
