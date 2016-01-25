#!/bin/bash

. /etc/profile

export SGE_ROOT=/var/lib/gridengine

PYTHON_COMMAND = ${PYTHONCOMMAND}

PIPELINE_ARGUMENTS=$*

#SHARED_FILESYSTEM_ROOT=/home/users/${USER}
SHARED_FILESYSTEM_ROOT=${HOME}

GENERAL_QUEUE=main.q

ASSEMBLER_QUEUE=large_host.q

RUNNAME=${RUNNAME}

LOGFILE=${SHARED_FILESYSTEM_ROOT}/data/pipeline_script${JOBNAME}.log

mkdir -p ${SHARED_FILESYSTEM_ROOT}/data/${RUNNAME}

echo ------------- NEW RUN ------------- >> $LOGFILE
date >> $LOGFILE
printenv >> $LOGFILE
printenv

${PYTHONCOMMAND} 2>&1 >> $LOGFILE 
