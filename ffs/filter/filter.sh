#!/bin/bash
# usage: ./filter.sh [OPTIONS]
# On can switch between NFS-DL or FFS computation with the following option:
#     nfsdl=[0|1]                   (default 0) FFS (0) or NFS-DL (1) 
# Mandatory options:
#     name=<name>
#     rels=/path/to/rels
#     cadobuild=/path/to/cado/bin
#     param=/path/to/param/file
# Mandatory options for NFS-DL (for FFS they are not used):
#     ell=NN                        group size (modulus for linear algebra)
#     smexp=NN                      exponent for the Shirokauer maps
# Other options:
#     wdir=path/to/output/directory (default ./<name>.filter.`date`)
#             if wdir already exists, continue previous computation
#     tidy=[0|1]                    (default 0)
#     covernmax=nn.nn               (default 100)
#     excess=nn                     (default 0 for FFS, deg(poly_alg) for NFS-DL)
#     req-excess=nn.nn              (default undefined)
#     maxlevel=nn                   (default 30)
#     addfullcol=[0|1]              (default 0)
#     badideals=file                (default "")
#     verbose=[0|1]                 (default 0)
# ex:
#  ./filter.sh name=ffs809 rels=/local/rsa768/ffs809/rels cadobuild=$HOME/cado-nfs/build/`hostname` param=/local/rsa768/ffs809/param.2.809
#
# TODO: addfullcol should be 1 by default for NFS-DL, 0 for FFS

check_error () {
  if [ $1 -ne 0 ] ; then
    echo "Error: $2"
    exit 1
  fi
}

run_cmd () {
  __cmd=$1
  __log1=$2
  __log2=$3
  __v=$4
  __donefile=$5

  __name=`basename ${__donefile} "_done"`

  echo "Start ${__name}"

  if [ "x${__log1}" = "x${__log2}" ] ; then
    if [ "x${__v}" = "x1" ] ; then
      echo "${__cmd} > ${__log1} 2>&1"
    fi

    ${__cmd} > ${__log1} 2>&1
    check_error "$?" "see ${__log1} for more info."
    echo "${__cmd} > ${__log1} 2>&1" >> ${CMDFILE}
  else
    if [ "x${__v}" = "x1" ] ; then
      echo "${__cmd} > ${__log1} 2> ${__log2}"
    fi

    ${__cmd} > ${__log1} 2>${__log2}
    check_error "$?" "see ${__log2} for more info."
    echo "${__cmd} > ${__log1} 2>${__log2}" >> ${CMDFILE}
  fi

  echo "End ${__name}"
  touch ${__donefile}
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
        declare RELDIR=`readlink -f $value`
      elif [ "x${p}" = "xcadobuild" ] ; then
        declare CADO_BUILD=`readlink -f $value`
      elif [ "x${p}" = "xparam" ] ; then
        declare ORIGINAL_PARAMFILE=`readlink -f $value`
      elif [ "x${p}" = "xwdir" ] ; then
        declare DIR=`readlink -f $value`
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
      elif [ "x${p}" = "xverbose" ] ; then
        declare VERBOSE="$value"
      elif [ "x${p}" = "xaddfullcol" ] ; then
        if [ "x$value" = "x1" ] ; then
          declare ADDFULLCOL="-addfullcol"
        fi
      elif [ "x${p}" = "xbadideals" ] ; then
        declare BADIDEALS="-badideals $value"
      elif [ "x${p}" = "xnfsdl" ] ; then
        if [ "x$value" = "x1" ] ; then
          declare NFSDL="1"
        fi
      elif [ "x${p}" = "xell" ] ; then
        declare MOD_ELL="$value"
      elif [ "x${p}" = "xsmexp" ] ; then
        declare SM_EXP="$value"
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
: ${EXCESS:="0"}
: ${REQ_EXCESS:="-1"}  # -1.0 means no parameter is passed in purge
: ${MAXLEVEL:="30"}
: ${BADIDEALS:=""}
: ${ADDFULLCOL:=""}
: ${VERBOSE:="0"}
: ${NFSDL:="0"}
: ${MOD_ELL:="<modulus>"}
: ${SM_EXP:="<sm_exp>"}


# Def variables to continue a filtering step that was stopped ######

DO_INIT="1"
DO_FREEREL="1"
DO_DUP1="1"
DO_DUP20="1"
DO_DUP21="1"
DO_PURGE="1"
DO_MERGE="1"
DO_REPLAY="1"
DO_SM="1"

# Init variables

mkdir -p ${DIR}

#For FFS case, recover the caracteristic from the param file
if [ "x${NFSDL}" != "x1" ] ; then
  GF=`grep "^gf=" ${ORIGINAL_PARAMFILE} | cut -d = -f 2`
else
  LPBA=`grep "^lpba" ${ORIGINAL_PARAMFILE} | cut -d " " -f 2`
  LPBR=`grep "^lpbr" ${ORIGINAL_PARAMFILE} | cut -d " " -f 2`
  if [ "x${EXCESS}" = "x0" ] ; then
    EXCESS=`grep -m 1 "^c[0-9]:" ${ORIGINAL_PARAMFILE} | cut -c 2`
  fi
fi

CMDFILE="${DIR}/${NAME}.cmd"

FILELIST_DUP1="${DIR}/${NAME}.filelist.dup1"
FILELIST_DUP2="${DIR}/${NAME}.filelist.dup2"
SUBDIRLIST="${DIR}/${NAME}.subdirlist"
NODUPDIR="${DIR}/${NAME}.nodup"
PARAMFILE="${ORIGINAL_PARAMFILE}"

LOGF="${DIR}/${NAME}.freerels.log"
LOGD1="${DIR}/${NAME}.dup1.log"
LOGD20="${DIR}/${NAME}.dup2_0.log"
LOGD21="${DIR}/${NAME}.dup2_1.log"
LOGP="${DIR}/${NAME}.purge.log"
LOGM="${DIR}/${NAME}.merge.log"
LOGR="${DIR}/${NAME}.replay.log"
LOGSM="${DIR}/${NAME}.sm.log"

NRELSFILE="${DIR}/${NAME}.nrels"

FREERELSFILE="${DIR}/${NAME}.freerels"
FREEGZFILE="${DIR}/${NAME}.freerels.gz"
RENUMBERFILE="${DIR}/${NAME}.renumber"
RELSFILE="${DIR}/${NAME}.rels.purged.gz"
DELRELSFILE="${DIR}/${NAME}.rels.deleted"
HISFILE="${DIR}/${NAME}.merge.his"
INDEXFILE="${DIR}/${NAME}.replay.index"
SMFILE="${DIR}/${NAME}.sm"

PREFIX_MATRIX="${DIR}/${NAME}.matrix"
INDEX_ID_MERGE="${DIR}/${NAME}.ideals"

#For FFS case, use dup-ffs-f*. For NFS-DL use dup2
if [ "x${NFSDL}" != "x1" ] ; then
BIN_FREE="${CADO_BUILD}/ffs/f${GF}/freerels";
BIN_DUP2="${CADO_BUILD}/filter/dup2-ffs-f${GF}";
else
BIN_FREE="${CADO_BUILD}/sieve/freerel";
BIN_DUP2="${CADO_BUILD}/filter/dup2";
fi
BIN_DUP1="${CADO_BUILD}/filter/dup1";
BIN_PURGE="${CADO_BUILD}/filter/purge";
BIN_MERGE="${CADO_BUILD}/filter/merge-dl";
BIN_REPLAY="${CADO_BUILD}/filter/replay-dl";
BIN_SM="${CADO_BUILD}/linalg/sm";
check_file_exists "${BIN_FREE}"
check_file_exists "${BIN_DUP1}"
check_file_exists "${BIN_DUP2}"
check_file_exists "${BIN_PURGE}"
check_file_exists "${BIN_MERGE}"
check_file_exists "${BIN_REPLAY}"
if [ "x${NFSDL}" = "x1" ] ; then
check_file_exists "${BIN_SM}"
fi

INITDONE="${DIR}/${NAME}.init_done"
FREEDONE="${DIR}/${NAME}.freerels_done"
DUP1DONE="${DIR}/${NAME}.dup1_done"
DUP20DONE="${DIR}/${NAME}.dup2_0_done"
DUP21DONE="${DIR}/${NAME}.dup2_1_done"
PURGEDONE="${DIR}/${NAME}.purge_done"
MERGEDONE="${DIR}/${NAME}.merge_done"
REPLAYDONE="${DIR}/${NAME}.replay_done"
SMDONE="${DIR}/${NAME}.sm_done"

if [ -e ${INITDONE} ] ; then
  DO_INIT="0"
  if [ -e ${FREEDONE} ] ; then
    DO_FREEREL="0"
    if [ -e ${DUP1DONE} ] ; then
      DO_DUP1="0"
      if [ -e ${DUP20DONE} ] ; then
        DO_DUP20="0"
        if [ -e ${DUP21DONE} ] ; then
          DO_DUP21="0"
          if [ -e ${PURGEDONE} ] ; then
            DO_PURGE="0"
            if [ -e ${MERGEDONE} ] ; then
              DO_MERGE="0"
              if [ -e ${REPLAYDONE} ] ; then
                DO_REPLAY="0"
                if [ -e ${SMDONE} ] ; then
                  DO_SM="0"
                fi
              fi
            fi
          fi
        fi
      fi
    fi
  fi
fi

###### INIT ######
if [ "${DO_INIT}" -eq "1" ] ; then

  echo "Start initialisation."
 
  ls -1 "${RELDIR}/" > ${FILELIST_DUP1} # all files in RELDIR goes to FILELIST
  check_error "$?" "could not construct list of relation files from ${RELDIR}."
  NB_RELFILE=`wc -l ${FILELIST_DUP1} | cut -d " " -f 1`
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

  touch ${INITDONE}
  echo "End initialisation."
else
  echo "initialisation already done."
fi
##################

###### FREERELS ######
if [ "${DO_FREEREL}" -eq "1" ] ; then

#For FFS case, use param option is different
if [ "x${NFSDL}" != "x1" ] ; then
  argsf1="${PARAMFILE} -renumber ${RENUMBERFILE} "
  argsf2=""
else
  argsf1="-poly ${PARAMFILE} -renumber ${RENUMBERFILE} "
  argsf2="-lpba ${LPBA} -lpbr ${LPBR} "
fi
  argsf3="-out ${FREEGZFILE} ${BADIDEALS} ${ADDFULLCOL} "
  CMD="${BIN_FREE} $argsf1 $argsf2 $argsf3"
  run_cmd "${CMD}" "${LOGF}" "${LOGF}" "${VERBOSE}" "${FREEDONE}"
else
  echo "freerels already done."
fi
######################


###### DUP1 ######
if [ "${DO_DUP1}" -eq "1" ] ; then
  argsd1="-n 1 -abhexa -filelist ${FILELIST_DUP1} -basepath ${RELDIR}"
  argsd2="-prefix ${NAME}.nodup.rels"

  CMD="${BIN_DUP1} $argsd1 $argsd2 -out ${NODUPDIR}"
  run_cmd "${CMD}" "${LOGD1}" "${LOGD1}" "${VERBOSE}" "${DUP1DONE}"
else
  echo "dup1 already done."
fi
##################
 
ls -1 "${NODUPDIR}/0/" > ${FILELIST_DUP2}
argsd2="-filelist ${FILELIST_DUP2} -poly ${PARAMFILE} -renumber ${RENUMBERFILE} "

###### DUP2_0 ######
if [ "${DO_DUP20}" -eq "1" ] ; then
  NRELSDUP20=`grep "^# slice 0" ${LOGD1} | cut -d " " -f 5`
  argsd20="-nrels ${NRELSDUP20} -basepath ${NODUPDIR}/0 "
  if [ "x${NFSDL}" != "x1" ] ; then
    CMD="${BIN_DUP2} $argsd2 $argsd20"
  else
    CMD="${BIN_DUP2} -dl $argsd2 $argsd20"
  fi
  run_cmd "${CMD}" "${LOGD20}" "${LOGD20}" "${VERBOSE}" "${DUP20DONE}"
else
  echo "dup2_0 already done."
fi
####################

###### DUP2_1 ######
if [ "${DO_DUP21}" -eq "1" ] ; then
  NRELSDUP21=`grep "^# slice 1" ${LOGD1} | cut -d " " -f 5`
  argsd21="-nrels ${NRELSDUP21} -basepath ${NODUPDIR}/1 "
  if [ "x${NFSDL}" != "x1" ] ; then
    CMD="${BIN_DUP2} $argsd2 $argsd21"
  else
    CMD="${BIN_DUP2} -dl $argsd2 $argsd21"
  fi
  run_cmd "${CMD}" "${LOGD21}" "${LOGD21}" "${VERBOSE}" "${DUP21DONE}"
else
  echo "dup2_1 already done."
fi
####################

NB0=`grep "^At the end: [0-9]" ${LOGD20} | cut -d " " -f 4`
NB1=`grep "^At the end: [0-9]" ${LOGD21} | cut -d " " -f 4`
NBREL=`expr ${NB0} + ${NB1}`
echo ${NBREL} > ${NRELSFILE}
echo "${NBREL} unique relations remaining."

NBPR=`grep nprimes ${LOGD20} | cut -d "=" -f 2`
#MIN=`grep min_index ${LOGD20} | cut -d "=" -f 2`
MIN="0" # FIXME MIN=ceil(2*minlim/log(minlim)) (in the integers case, not true
#in the polynomials case (as in FFS))

###### PURGE ######
if [ "${DO_PURGE}" -eq "1" ] ; then
  argp0="-out ${RELSFILE} -basepath ${NODUPDIR} -subdirlist ${SUBDIRLIST} "
  argp1="-filelist ${FILELIST_DUP2} -keep ${EXCESS} "
  argp2="-outdel ${DELRELSFILE} "
  argp3="-nrels ${NBREL} -nprimes ${NBPR} -minindex ${MIN} "
  if [ "x${REQ_EXCESS}" = "x-1" ] ; then
    argp4=""
  else
    argp4="-required_excess ${REQ_EXCESS} "
  fi

  CMD="${BIN_PURGE} $argp0 $argp1 $argp2 $argp3 $argp4"
  run_cmd "${CMD}" "${LOGP}" "${LOGP}" "${VERBOSE}" "${PURGEDONE}"
else
  echo "purge already done."
fi
###################

STAT_END_PURGE=`grep "^nrels" ${LOGP}`
echo "${STAT_END_PURGE} ."

###### MERGE ######
if [ "${DO_MERGE}" -eq "1" ] ; then
  argm0="-out ${HISFILE} -mat ${RELSFILE} -forbw 3 -coverNmax ${COVERNMAX} "
  argm1="-keep ${EXCESS} -maxlevel ${MAXLEVEL} -skip 0 "

  CMD="${BIN_MERGE} $argm0 $argm1"
  run_cmd "${CMD}" "${LOGM}" "${LOGM}" "${VERBOSE}" "${MERGEDONE}"

  RELS_AFTER_MERGE=`grep "Final mat" ${LOGM} | cut -d " " -f 4 | cut -d = -f 2`
  echo "Matrix has ${RELS_AFTER_MERGE} rows."
else
  echo "merge already done."
fi
###################

###### REPLAY ######
if [ "${DO_REPLAY}" -eq "1" ] ; then
  if [ "x${NFSDL}" != "x1" ] ; then
    argsr0="--noindex -purged ${RELSFILE} -his ${HISFILE} "
  else
    argsr0="-index ${INDEXFILE} -purged ${RELSFILE} -his ${HISFILE} "
  fi
  argsr1="-out ${PREFIX_MATRIX} -ideals ${INDEX_ID_MERGE} -skip 0"

  CMD="${BIN_REPLAY} $argsr0 $argsr1 $argsr2"
  run_cmd "${CMD}" "${LOGR}" "${LOGR}" "${VERBOSE}" "${REPLAYDONE}"
else
  echo "replay already done."
fi
####################


######## SM ########   # only for NFS-DL
if [ "x${NFSDL}" = "x1" ] ; then
  if [ "${DO_SM}" -eq "1" ] ; then
    echo "${MOD_ELL}" > ${DIR}/${NAME}.q
    argssm0="-poly ${PARAMFILE} -purged ${RELSFILE} -index ${INDEXFILE} "
    argssm1="-out ${SMFILE} -gorder ${MOD_ELL} -smexp ${SM_EXP}"
    CMD="${BIN_SM} $argssm0 $argssm1 "
    run_cmd "${CMD}" "${LOGSM}" "${LOGSM}" "${VERBOSE}" "${SMDONE}"
  else
    echo "sm already done."
  fi
fi
####################


# Help: output the command line needed to recontruct all logarithms from the
# ones computed by linear algebra

BIN_RECONSTRUCT="${CADO_BUILD}/filter/reconstructlog-ffs-f${GF}";
OUTLOG="${DIR}/${NAME}.logarithms.values"
argsre0="-ideals ${INDEX_ID_MERGE} -relsdel ${DELRELSFILE} -nrels ${NBREL}"
argsre1="-relspurged ${RELSFILE} -renumber ${RENUMBERFILE} -poly ${PARAMFILE}"
argsre2="-out ${OUTLOG} -log <file> -q ${MOD_ELL}"
echo "${BIN_RECONSTRUCT} $argsre0 $argsre1 $argsre2 "

###### If tidy is asked ######

if [ "${TIDY}" -eq "1" ] ; then
  rm -r ${NODUPDIR} ${INVERTDIR}
  rm ${PARAMFILE} ${FILELIST} ${SUBDIRLIST}
  rm "${RELSFILE}" "${NRELSFILE}"
  rm "${DELRELSFILE}"
  rm "${HISFILE}"
  rm "${INDEX_ID_MERGE}"
fi
