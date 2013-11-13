# Run inside the parameter directory

PARAMDIR=.
if [ -n "$1" ]
then
  PARAMDIR="$1"
fi

if ! test -d "$PARAMDIR" || ! test -f "$PARAMDIR/params.c91"
then
  echo Please specify the parameter directory
  exit 1
fi

for S in degree qrange I polsel_P polsel_maxnorm polsel_admax polsel_adrange rlim alim lpbr lpba
do
  OUTPUTFILE="$S.png"
  echo "Creating $OUTPUTFILE"
  grep "^$S=" "$PARAMDIR"/params.c[1-9][0-9] "$PARAMDIR"/params.c1[0-9][0-9] | sed "s/^.*params.c//; s/:$S=/ /" | cut -d " " -f 1,2 >| /tmp/$S.data
  echo "set terminal png; set log y 2; set output \"$OUTPUTFILE\"; plot \"/tmp/$S.data\" " | gnuplot
done
echo "Creating montage.png"
montage -geometry 640x480 *.png montage.png
