#!/bin/bash

#Current Version is for "Direct" only

if [ $# -eq 0 ]
then
  echo >&2 "ERROR: $0 need at least 2 options. Try \"$0 --help\" for help"
  exit 0
fi

if [ $# -eq 1 ] 
then
  if [ $1 != '--help' ]
  then
    echo >&2 "ERROR: $0 need at least 2 options. Try \"$0 --help\" for help"
    exit 0
  else
    echo >&2 "Usage: $0 POSFILE_TO_READ VARIABLE_TO_READ [OPTION1 [OPTION2]]"
    echo >&2 "- VARAIBLE_TO_READ="
    echo >&2 "- lc: Lattice Constant"
    echo >&2 "- lv: Lattice Vector"
    echo >&2 "- sc: Species Count"
    echo >&2 "- sn: Species Name"
    echo >&2 "- sa/aps: Atoms per Species"
    echo >&2 "- ac: Atom Count"
    echo >&2 "- an: Atom Name"
    echo >&2 "- ap: Atom Position"
    echo >&2 "For Detailed Help, try \"$0 VARIABLE_TO_READ\"."
    exit 0
  fi
else
  if [ ! -f $1 ]
  then
    echo >&2 "ERROR: POSFILE file doesn't exist!"
    exit 1
  fi
  poscar_file=$1
  shift
  case $1 in 
  lattconst|lc)
    sed -n "2p" ${poscar_file} | sed 's/ //g'  
    ;;
  lattvec|lv)
    #relative lattice vector, with a scaling of lattconst
    #must specify option1 & option2
    if [ $# -ne 3 ]
    then 
      echo >&2 "ERROR: Need OPTION1 (Latt. Vec. No.) & OPTION2 (which axis: x~z=1~3)"
      echo >&2 "OPTION1 or OPTION2 could be taken as 0=all"
      echo >&2 "if OPTION1 == all, OPTION2 is automatically set to all,"
      echo >&2 " and all lattice vectors will be print."
      echo >&2 "if OPTION2 == all, the entire lattice vector will be print (all 3 axis)."
      exit 1
    fi
    isvalid='^[0-3]$'
    if ! [[ $2 =~ $isvalid ]]
    then
      echo >&2 "ERROR: OPTION1 must be a number between 0~3."
      exit 1
    fi
    if ! [[ $3 =~ $isvalid ]]
    then
      echo >&2 "ERROR: OPTION2 must be a number between 0~3."
      exit 1
    fi
    if [ $2 -eq 0 ]
    then
      sed -n "3,5p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g'
    elif [ $3 -eq 0 ]
    then
      sed -n "$((2+$2))p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g'
    else
      sed -n "$((2+$2))p" ${poscar_file} | awk "{print $"$3";}"
    fi
    ;; 
  speciescount|sc)
    sed -n "6p" ${poscar_file} | awk "{print NF;}"    
    ;;
  speciesname|sn)
    #must specify option1
    if [ $# -ne 2 ]
    then
      echo >&2 "ERROR: Need OPTION1 (Species No. beginning with 1; 0 = all)"
      exit 1
    fi
    scount=`sed -n "6p" ${poscar_file} | awk "{print NF;}"`
    sno=$2
    isnumber='^[0-9]+$'
    if ! [[ $sno =~ $isnumber ]]
    then
      echo >&2 "ERROR: OPTION1 (Species No.) need to be a number!"
      exit 1
    fi
    if [ $sno -gt $scount ]
    then
      echo >&2 "ERROR: OPTION1 (Species No.) = $sno > $scount = Species Count."
      exit 1
    fi
    if [ $sno -eq 0 ]
    then
      sed -n "6p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g'
    else
      sed -n "6p" ${poscar_file} | awk "{print $"$sno";}"
    fi
    ;;
  atomcount|ac)
    sed -n "7p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g' | sed 's/ \+/+/g' | bc
    ;;
  atomname|an)
    #must specify option1
    if [ $# -ne 2 ]
    then
      echo >&2 "ERROR: Need OPTION1 (Atom No. beginning with 1; 0 = all)"
      exit 1
    fi
    scount=`sed -n "7p" ${poscar_file} | awk "{print NF;}"`
    acount=`sed -n "7p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g' | sed 's/ \+/+/g' | bc`
    ano=$2
    isnumber='^[0-9]+$'
    if ! [[ $ano =~ $isnumber ]]
    then
      echo >&2 "ERROR: OPTION1 (Atom No.) need to be a number!"
      exit 1
    fi
    if [ $ano -gt $acount ]
    then
      echo >&2 "ERROR: OPTION1 (Atom No.) = $ano > $acount = Atom Count."
      exit 1
    fi

    #echo $scount $acount $ano
    # Find the species corresponding to atom #iano, for iano in 1-acount
    aupperlim=0
    scur=0
    for iano in `seq 1 $acount`
    do
      #echo $iano
      if [ $iano -gt $aupperlim ]
      then
        scur=$((scur+1))
        #echo $scur
        sacount=`sed -n "7p" ${poscar_file} | awk "{print $"$scur";}"`
        #echo $sacount
        aupperlim=$(($aupperlim+$sacount))
      fi
      if [ $ano -eq 0 ] || [ $iano -eq $ano ]
      then
        sed -n "6p" ${poscar_file} | awk "{printf $"$scur"; printf \" \";}"
      fi
    done
    echo

    ;;
  atomperspecies|aps|speciesatoms|sa)
    #must specify option1
    if [ $# -ne 2 ]
    then
      echo >&2 "ERROR: Need OPTION1 (Species No. beginning with 1; 0 = all)"
      exit 1
    fi
    scount=`sed -n "7p" ${poscar_file} | awk "{print NF;}"`
    sno=$2
    isnumber='^[0-9]+$'
    if ! [[ $sno =~ $isnumber ]]
    then
      echo >&2 "ERROR: OPTION1 (Species No.) need to be a number!"
      exit 1
    fi
    if [ $sno -gt $scount ]
    then
      echo >&2 "ERROR: OPTION1 (Species No.) = $sno > $scount = Species Count."
      exit 1
    fi
    if [ $sno -eq 0 ]
    then
      sed -n "7p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g'
    else
      sed -n "7p" ${poscar_file} | awk "{print $"$sno";}"
    fi
    ;;
  atomposition|ap)
    #must specify option1 & option 2
    if [ $# -ne 3 ]
    then
      echo >&2 "ERROR: Need OPTION1 (Atom. No.) & OPTION2 (which axis: x~z=1~3)"
      echo >&2 "OPTION1 or OPTION2 could be taken as 0=all"
      echo >&2 "if OPTION1 == all, OPTION2 is automatically set to all,"
      echo >&2 " and all atom positions will be print."
      echo >&2 "if OPTION2 == all, the entire atom position vector will be print (all 3 axis)."
      exit 1
    fi
    acount=`sed -n "7p" ${poscar_file} | sed 's/^ *\(.*\) *$/\1/g' | sed 's/ \+/+/g' | bc`
    ano=$2
    isnumber='^[0-9]+$'
    if ! [[ $ano =~ $isnumber ]]
    then
      echo >&2 "ERROR: OPTION1 (Atom No.) need to be a number!"
      exit 1
    fi
    if [ $ano -gt $acount ]
    then
      echo >&2 "ERROR: OPTION1 (Atom No.) = $ano > $acount = Atom Count."
      exit 1
    fi
    isvalid='^[0-3]$'
    if ! [[ $3 =~ $isvalid ]]
    then
      echo >&2 "ERROR: OPTION2 must be a number between 0~3."
      exit 1
    fi
    if [ $ano -eq 0 ]
    then
      sed -n "9,$((8+$acount))p" ${poscar_file} | awk '{printf "%20.16f %20.16f %20.16f\n", $1, $2, $3;}'
    elif [ $3 -eq 0 ]
    then
      sed -n "$((8+$ano))p" ${poscar_file} | awk '{printf "%20.16f %20.16f %20.16f\n", $1, $2, $3;}'
    else
      sed -n "$((8+$ano))p" ${poscar_file} | awk "{print $"$3";}"
    fi
    ;;
  *)
    echo >&2 "ERROR: Cannot identify VARIABLE_TO_READ: $1"
    exit 1
    ;;
  esac

fi
