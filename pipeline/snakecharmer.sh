#!/usr/bin/env bash
#BSUB -J squint
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 

set -o nounset -o pipefail -o errexit -x

args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R
"{params.memory} span[hosts=1] " -n {threads} '
    
## make sure that following are all in your path

# salmon
# gffread

#### execute snakemake ####

snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 24 \
    --resources all_threads=24 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config.yaml 
