#!/bin/bash
# usage: ./filter.sh [OPTIONS]
# Mandatory options:
#     name=<name> 
#     rels=/path/to/rels 
#     cadobuild=/path/to/cado/bin 
#     param=/path/to/param/file
# Other options:
#     wdir=path/to/output/directory (default ./<name>.filter.`date`)
#             if wdir already exists, continue previous computation
#     tidy=[0|1]                    (default 0)
#     covernmax=nn.nn               (default 100)
#     excess=nn                     (default 1)
#     req-excess=nn.nn              (default undefined)
#     maxlevel=nn                   (default 30)
# ex: 
#  ./filter.sh name=ffs809 rels=/local/rsa768/ffs809/rels cadobuild=$HOME/cado-nfs/build/`hostname` param=/local/rsa768/ffs809/param.2.809

check_error () {
  if [ $1 -ne 0 ] ; then
    echo "Error: $2"
    exit 1
  fi
}

check_file_exists () {
  if [ ! -f "$1" ] ; then
    echo "$1 does not exist or is not a regular file."
    exit 1
  fi
}

check_dir_exists () {
  if [ ! -d "$1" ] ; then
    echo "$1 does not exist or is not a directory."
    exit 1
  fi
}

for i in "$@" ; do 
  case "$i" in
    *=*)
      shift
      value=${i#*=}
      p=${i%%=*}
      if [ "x${p}" = "xname" ] ; then
        declare NAME="$value"
      elif [ "x${p}" = "xrels" ] ; then
        declare RELDIR="$value"
      elif [ "x${p}" = "xcadobuild" ] ; then
        declare CADO_BUILD="$value"
      elif [ "x${p}" = "xparam" ] ; then
        declare ORIGINAL_PARAMFILE="$value"
      elif [ "x${p}" = "xwdir" ] ; then
        declare DIR="$value"
      elif [ "x${p}" = "xtidy" ] ; then
        declare TIDY="$value"
      elif [ "x${p}" = "xcovernmax" ] ; then
        declare COVERNMAX="$value"
      elif [ "x${p}" = "xexcess" ] ; then
        declare EXCESS="$value"
      elif [ "x${p}" = "xreq-excess" ] ; then
        declare REQ_EXCESS="$value"
      elif [ "x${p}" = "xmaxlevel" ] ; then
        declare MAXLEVEL="$value"
      else
        echo "Invalid option: $i"
        exit 1
      fi
      ;;
    *)
      shift
      echo "Invalid option: $i"
      exit 1;
  esac
done

# check that madatory options are set and non null
: ${NAME:?Option name=<name> is mandatory}
: ${RELDIR:?Option rels=/path/to/rels is mandatory}
: ${CADO_BUILD:?Option cadobuild=/path/to/cado/bin is mandatory}
: ${ORIGINAL_PARAMFILE:?Option param=/path/to/param/file is mandatory}

check_dir_exists "${RELDIR}"
check_file_exists "${ORIGINAL_PARAMFILE}"

# put default values of options that are non set or null
: ${DIR:="${NAME}.filter.`date +%y%m%d-%H%M`"}
: ${TIDY:="0"}
: ${COVERNMAX:="100"}
: ${EXCESS:="1"} 
: ${REQ_EXCESS:="-1"}  # -1.0 means no parameter is passed in purge
: ${MAXLEVEL:="30"}


# Def variables to continue a filtering step that was stopped ######

DO_INIT="1"
DO_DUP="1"
DO_INVERT="1"
DO_PURGE="1"
DO_MERGE="1"
DO_REPLAY="1"

# Init variables

mkdir -p ${DIR}

GF=`grep "^gf=" ${ORIGINAL_PARAMFILE} | cut -d = -f 2`

#invert param file and rel if pol0 is the algebraic side
RAT_SIDE=`grep -c "pol0=[0-9a-zA-Z]*,[0-9a-zA-Z]*," ${ORIGINAL_PARAMFILE}`
if [ ${RAT_SIDE} -eq "0" ] ; then 
  DO_INVERT="0"
fi

FILELIST="${DIR}/${NAME}.filelist"
SUBDIRLIST="${DIR}/${NAME}.subdirlist"
NODUPDIR="${DIR}/${NAME}.nodup"

if [ "${DO_INVERT}" -eq "1" ] ; then
  INVERTDIR="${DIR}/${NAME}.invert"
  PARAMFILE="${DIR}/${NAME}.param.inv"
else
  INVERTDIR="${NODUPDIR}"
  PARAMFILE="${ORIGINAL_PARAMFILE}"
fi

LOGD1="${DIR}/${NAME}.dup1.log"
LOGD20="${DIR}/${NAME}.dup2-0.log"
LOGD21="${DIR}/${NAME}.dup2-1.log"
LOGP="${DIR}/${NAME}.purge.log"
LOGM="${DIR}/${NAME}.merge.log"
LOGR="${DIR}/${NAME}.replay.log"

NRELSFILE="${DIR}/${NAME}.nrels"

RELSFILE="${DIR}/${NAME}.rels.purged.gz"
INDEX_ID_PURGE="${DIR}/${NAME}.ideals.tmp"
DELRELSFILE="${DIR}/${NAME}.rels.deleted"
HISFILE="${DIR}/${NAME}.merge.his"

PREFIX_MATRIX="${DIR}/${NAME}.matrix"
MERGEDRELSFILE="${DIR}/${NAME}.rels.merged"
INDEX_ID_MERGE="${DIR}/${NAME}.id.replay.index"

BIN_DUP1="${CADO_BUILD}/filter/dup1";
BIN_DUP2="${CADO_BUILD}/filter/dup2";
BIN_PURGE="${CADO_BUILD}/filter/purge-ffs-f${GF}";
BIN_MERGE="${CADO_BUILD}/filter/merge-ffs";
BIN_REPLAY="${CADO_BUILD}/filter/replay-ffs";
check_file_exists "${BIN_DUP1}"
check_file_exists "${BIN_DUP2}"
check_file_exists "${BIN_PURGE}"
check_file_exists "${BIN_MERGE}"
check_file_exists "${BIN_REPLAY}"

INITDONE="${DIR}/${NAME}.init_done"
DUPDONE="${DIR}/${NAME}.dup_done"
if [ "${DO_INVERT}" -eq "1" ] ; then
  INVERTDONE="${DIR}/${NAME}.invert_done"
else
  INVERTDONE=""
fi
PURGEDONE="${DIR}/${NAME}.purge_done"
MERGEDONE="${DIR}/${NAME}.merge_done"
REPLAYDONE="${DIR}/${NAME}.replay_done"

if [ -e ${INITDONE} ] ; then 
  DO_INIT="0"
  if [ -e ${DUPDONE} ] ; then 
    DO_DUP="0"
    if [ -e ${INVERTDONE} ] ; then 
      DO_INVERT="0"
      if [ -e ${PURGEDONE} ] ; then 
        DO_PURGE="0"
        if [ -e ${MERGEDONE} ] ; then 
          DO_MERGE="0"
          if [ -e ${REPLAYDONE} ] ; then 
            DO_REPLAY="0"
          fi
        fi
      fi
    fi
  fi
fi

###### INIT ######
if [ "${DO_INIT}" -eq "1" ] ; then

  echo "Start initialisation."
  
  ls -1 "${RELDIR}/" > ${FILELIST} # all files in RELDIR goes to FILELIST
  check_error "$?" "could not construct list of relation files from ${RELDIR}."
  NB_RELFILE=`wc -l ${FILELIST} | cut -d " " -f 1`
  if [ ${NB_RELFILE} -eq 0 ] ; then
    echo "The directory ${RELDIR} contains no relations file."
    exit 1
  fi
  echo "  create filelist with ${NB_RELFILE} relations files."

  echo "0" > ${SUBDIRLIST}      # SUBDIRLIST is just 2 lines: 0 and 1
  echo "1" >> ${SUBDIRLIST}
  echo "  create subdirlist."

  mkdir -p ${NODUPDIR}/0 ${NODUPDIR}/1
  check_error "$?" "could not create ${NODUPDIR}/0 and ${NODUPDIR}/1."

  # Transform param file from alg=0 rat=1 to rat=0 alg=1 for purge
  #TODO invert I & J ? sqside ? S ?
if [ "${DO_INVERT}" -eq "1" ] ; then
  mkdir -p ${INVERTDIR}/0 ${INVERTDIR}/1
  check_error "$?" "could not create ${INVERTDIR}/0 and ${INVERTDIR}/1."

  sed -e 's/0=/42=/gI' -e 's/1=/0=/gI' ${ORIGINAL_PARAMFILE} > ${DIR}/${NAME}.tmp
  sed -e 's/42=/1=/gI' ${DIR}/${NAME}.tmp > ${PARAMFILE}
  rm ${DIR}/${NAME}.tmp
  echo "  invert alg and rat side of param file."
fi

  touch ${INITDONE}
  echo "End initialisation."
else
  echo "initialisation already done."
fi

###### DUP ######
if [ "${DO_DUP}" -eq "1" ] ; then

  echo "Start duplicates."

  argsd1="-n 1 -abhexa -filelist ${FILELIST} -basepath ${RELDIR} "

  echo -n "  dup1..."
  ${BIN_DUP1} $argsd1 -out ${NODUPDIR} > ${LOGD1} 2>&1
  check_error "$?" "see ${LOGD1} for more info."
  echo "done"
  
  ESTIMATED_NB_REL="`tail -n 1 ${LOGD1} | cut -d \" \" -f 3`"
  argsd2="-abhexa -K ${ESTIMATED_NB_REL} -filelist ${FILELIST} "

  echo -n "  dup2..."
  ${BIN_DUP2} $argsd2 -out ${NODUPDIR}/0 -basepath ${NODUPDIR}/0 > ${LOGD20} 2>&1
  check_error "$?" "see ${LOGD20} for more info."
  ${BIN_DUP2} $argsd2 -out ${NODUPDIR}/1 -basepath ${NODUPDIR}/1 > ${LOGD21} 2>&1
  check_error "$?" "see ${LOGD21} for more info."
  echo "done"

  NB0=`tail -n 3 ${LOGD20} | head -n 1 | cut -d " " -f 6`
  NB1=`tail -n 3 ${LOGD21} | head -n 1 | cut -d " " -f 6`
  NBREL=`expr ${NB0} + ${NB1}`
  echo ${NBREL} > ${NRELSFILE}

  touch ${DUPDONE}
  echo "End duplicates. ${NBREL} unique relations remaining."
else
  echo "duplicates already done."
fi

###### invert ######
if [ "${DO_INVERT}" -eq "1" ] ; then

  echo "Starting to invert relations files."

  for file in `cat ${FILELIST}` ; do
    gzip -dc ${NODUPDIR}/0/$file | awk -f invert_alg_rat.awk | gzip -c --fast > ${INVERTDIR}/0/$file
    gzip -dc ${NODUPDIR}/1/$file | awk -f invert_alg_rat.awk | gzip -c --fast > ${INVERTDIR}/1/$file
  done

  touch ${INVERTDONE}
  echo "End invert relations files."
else
  echo "invert relations files already done."
fi

###### PURGE ######

if [ "${DO_PURGE}" -eq "1" ] ; then
  echo "Start purge."

  argp0="-poly ${PARAMFILE} -out ${RELSFILE} -basepath ${INVERTDIR} "
  argp1="-subdirlist ${SUBDIRLIST} -filelist ${FILELIST} -keep ${EXCESS} "
  argp2="-sos ${INDEX_ID_PURGE} -outdel ${DELRELSFILE} "
  if [ "x${REQ_EXCESS}" = "x-1" ] ; then
    argp3=""
  else
    argp3="-required_excess ${REQ_EXCESS} "
  fi
  NBREL=`cat ${NRELSFILE}`

  ${BIN_PURGE} $argp0 -nrels ${NBREL} $argp1 $argp2 $argp3 > $LOGP 2>&1
  check_error "$?" "see ${LOGP} for more info."

  touch ${PURGEDONE}

  RELS_AFTER_PURGE=`grep NROWS ${LOGP} | cut -d " " -f 1 | cut -d : -f 2`
  echo "End purge. ${RELS_AFTER_PURGE} relations remaining."
else
  echo "purge already done."
fi

###### MERGE ######

if [ "${DO_MERGE}" -eq "1" ] ; then
  echo "Start merge."

  argm0="-out ${HISFILE} -mat ${RELSFILE} -forbw 3 -coverNmax ${COVERNMAX} "
  argm1="-keep ${EXCESS} -maxlevel ${MAXLEVEL} -skip 0 " 

  ${BIN_MERGE} $argm0 $argm1 > $LOGM 2>&1
  check_error "$?" "see ${LOGM} for more info."

  touch ${MERGEDONE}

  RELS_AFTER_MERGE=`grep "Final mat" ${LOGM} | cut -d " " -f 4 | cut -d = -f 2`
  echo "End merge. Matrix has ${RELS_AFTER_MERGE} rows."
else
  echo "merge already done."
fi

###### REPLAY ######

if [ "${DO_REPLAY}" -eq "1" ] ; then

  echo "Start replay."

  argsr0="--noindex -purged ${RELSFILE} -his ${HISFILE} "
  argsr1="-out ${PREFIX_MATRIX} -ideals ${INDEX_ID_MERGE} "
  argsr2="-outdel ${MERGEDRELSFILE}"

  ${BIN_REPLAY} $argsr0 $argsr1 $argsr2 > $LOGR 2>&1
  check_error "$?" "see ${LOGR} for more info."

  touch ${REPLAYDONE}
  echo "End replay."
else
  echo "replay already done."
fi

###### If tidy is asked ######

if [ "${TIDY}" -eq "1" ] ; then
  rm -r ${NODUPDIR} ${INVERTDIR}
  rm ${PARAMFILE} ${FILELIST} ${SUBDIRLIST}
  rm "${RELSFILE}" "${NRELSFILE}"
  rm "${INDEX_ID_PURGE}" "${DELRELSFILE}"
  rm "${HISFILE}" "${MERGEDRELSFILE}"
  rm "${INDEX_ID_MERGE}"
fi
