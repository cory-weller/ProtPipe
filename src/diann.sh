#!/usr/bin/env bash


#### CHECK SINGULARITY #############################################################################
if ! command -v singularity &> /dev/null; then
    echo "INFO: singularity command not found, looking for singularity module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        if $(module avail singularity/3 2>&1 >/dev/null | grep -q 'No module'); then
            echo 'ERROR: singularity cannot be found. Recheck installation?'
            exit 1
        else
            echo 'INFO: module singularity found'
            module load singularity/3
        fi
    fi
else
    echo 'INFO: singularity command found'
fi

#### CHECK PYTHON3 #################################################################################
if ! command -v python3 &> /dev/null; then
    if command -v module &> /dev/null; then
        if $(module avail python/3 2>&1 >/dev/null | grep -q 'No module'); then
            echo "ERROR: python3 module does not exist. Recheck installation?"
            exit 1
        else
            module load python/3
        fi
    fi
else
    echo "INFO: python3 command found"
fi


#### RUN DIA-NN ####################################################################################

if grep -q -- '--skip' <<< $(echo $@); then
    echo 'INFO: Skipping DIA-NN due to --skip'
    # Remove --skip from $@
    for arg do
        shift
        [ "$arg" = "--skip" ] && continue
        set -- "$@" "$arg"
    done
else
    python3 src/diann.py $@
fi

