#!/bin/bash

# 3 versions for this system: normal, negative-y and positive-y version.
inpos=$1  #POSCAR file
posx=$(awk 'NR==9{print $1;}' ${inpos})
posy=$(awk 'NR==9{print $2;}' ${inpos})
posz=$(awk 'NR==9{print $3;}' ${inpos})
head -n 8 ${inpos}
awk -v "px=${posx}" -v "py=${posy}" -v "pz=${posz}"  'NR>8{
  ax=($1-px+3.5)-int($1-px+3.5)-0.5;
  ay=($2-py+3.5)-int($2-py+3.5)-0.5;
  az=($3-pz+3.5)-int($3-pz+3.5)-0.5;
  printf "%20.16f %20.16f %20.16f  %s\n", ax, ay, az, $4;
}' ${inpos}
