# Run inside the parameter directory

for S in degree qrange I polsel_P polsel_maxnorm polsel_admax polsel_adrange rlim alim lpbr lpba
do
  grep "^$S=" params.c[1-9][0-9] params.c[1-9][0-9][0-9] | sed "s/params.c//; s/:$S=/ /" | cut -d " " -f 1,2 >| /tmp/$S.data
  echo "set terminal png; set log y 2; set output \"$S.png\"; plot \"/tmp/$S.data\" " | gnuplot
done
montage -geometry 640x480 *.png montage.png
