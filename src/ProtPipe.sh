#!/usr/bin/env bash

singularity_vers='3'



#### CHECK SINGULARITY #############################################################################
if ! command -v singularity &> /dev/null; then
    echo "INFO: singularity command not found, looking for singularity module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        if $(module avail singularity/${singularity_vers} 2>&1 >/dev/null | grep -q 'No module'); then
            echo 'ERROR: singularity cannot be found. Recheck installation?'
            exit 1
        else
            echo 'INFO: module singularity found'
            module load singularity/${singularity_vers}
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
    python3 src/dia-nn.py $@
fi


#### R ANALYSIS ####################################################################################
# Check singularity image integrity
r_version='r/4.0:1.3'
r_sif="src/R.sif"
r_sif_md5_desired='b052311c61f29e1815d5b7436a35b8a2'
if [ ! -f "${r_sif}" ]; then
    echo "INFO: Pulling ${r_sif} from remote library://wellerca/${r_version}" 
    singularity pull ${r_sif} library://wellerca/${r_version}
else
    echo "INFO: ${r_sif} already exists, skipping download"
fi

r_sif_md5_actual=$(md5sum "${r_sif}" | awk '{print $1}')

if [ ! "${r_sif_md5_actual}" == "${r_sif_md5_desired}" ]; then
    echo "ERROR: ${r_sif} md5 sum does not pass check. Possibly corrupted? Delete and try again."
    exit 1
else
    echo "INFO: ${r_sif} md5 sum passes check"
fi

diann_r_file='src/counts_processing.R'
folder='output'
singularity exec --cleanenv -H ${PWD} ${r_sif} Rscript ${diann_r_file} \
--pgfile        ${folder}/report.pg_matrix.tsv \
--out           ${folder} \
--design        example/design_matrix.csv


