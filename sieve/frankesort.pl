#!/usr/bin/perl

while (<>)
{
  split;
  $T = shift(@_);
  if ($T eq "X" || $T eq "Y") 
  {
    print("$T ", join(" ",sort{hex($a) <=> hex($b)} @_), "\n");
  } else {
    print ("$T ", join(" ", @_), "\n");
  }
}
