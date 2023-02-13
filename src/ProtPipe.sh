#!/usr/bin/env bash

singularity_vers='3'

# Check singularity exists
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

# Check python3 exists
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

python3 src/dia-nn.py $@

# add singularity pull for R .sif

diann_r_file='src/diann_processing.R'
folder='exampleOutput'
singularity exec --cleanenv -H ${PWD} src/r-4.0.sif Rscript ${diann_r_file} \
--pro_input         ${folder}/report.pg_matrix.tsv \
--pep_input         ${folder}/report.gg_matrix.tsv \
--prefix            test \
--outdir            ${folder} \
--design_matrix     example/design_matrix.csv