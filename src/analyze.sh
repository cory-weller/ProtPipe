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



#### R ANALYSIS ####################################################################################
r_sif="src/R.sif"
r_sif_md5_desired='5360384b914c94c97b0aa39eb334e7f0'
if [ ! -f "${r_sif}" ]; then
    echo "INFO: dowloading ${r_sif} from Onedrive"
    wget -O "${r_sif}" 'https://onedrive.live.com/download?resid=77DD71E598E5B51B%2125156&authkey=!ACYXTqwaZOfpvAQ'
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

r_processing_script='src/counts_processing.R'
singularity exec --cleanenv -H ${PWD} ${r_sif} Rscript ${r_processing_script} $@
