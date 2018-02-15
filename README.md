
## Salmon QUantification of Introns in Nascent Transcripts (SQUINT)

Squint uses `salmon` to calculate premRNA and mature mRNA expression estimates. 

## Install

```bash
git clone git@github.com:kriemo/squint
```

## Usage

Squint is a `snakemake` pipeline that is written for use on `tesla`. To use `squint` 
perform the following steps.

1) `cd` to the `squint` directory and edit the `config.yaml` file to point the transcripts and genome to the correct 
   genome assemblies and transcriptomes. 

i.e.
```bash
TRANSCRIPTS:
   - "path/to/your/gtf.gtf"

GENOME:
   - "path/to/your/genome.fa"
```

2) Make sure you have `salmon` in your path

i.e. this works
```bash
salmon -h
```

if you do not have the salmon binaries you can either download them yourself from [here]((https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz)
make sure to add the `bin` folder to your PATH in your bashrc. 

3) move your `fastq.gz` files into the `data/raw_data` directory. `Snakemake` will autodetect these files and process
them, as long as the filenames end with `.fastq.gz`. If you need to remove adapters or perform any trimming, do this
before moving your fastq files into the `data/raw_data` folder.

4) Use snakemake to run all of your samples. `Snakemake` will take care of submitting all of your jobs and will track if
a file had successfully generated. First test that snakemake recognizes your files.

change working directory to the `pipeline` folder
```bash
cd pipeline
```

execute a dry snakemake run to see all of the commands and files that will be processed

```bash
snakemake -npr 
```

you should see many lines of output that show all of the commands that snakemake will run


```bash
processing the following libraries
library test
rule primary_transcript_gtf:
    input: /vol3/home/riemondy/Projects/shared_dbases/annotation/drosophila/Drosophila_melanogaster.BDGP6.84.gtf
    output: ../dbases/primary_transcripts.fa, ../dbases/primary_transcripts.gtf
    log: ../dbases/logs/generate_primary_transcript_db.txt
    reason: Missing output files: ../dbases/primary_transcripts.fa, ../dbases/primary_transcripts.gtf
    resources: all_threads=1

      python3 ../src/add_primary_transcripts.py         -i /vol3/home/riemondy/Projects/shared_dbases/annotation/drosophila/Drosophila_melanogaster.BDGP6.84.gtf         -r "gene"         -a "gene_id" > ../dbases/primary_transcripts.gtf

      gffread -w ../dbases/primary_transcripts.fa -g /vol3/home/riemondy/Projects/shared_dbases/genomes/drosophila/Drosophila_melanogaster.BDGP6.dna.toplevel.fa ../dbases/primary_transcripts.gtf

rule salmon_idx:
    input: ../dbases/primary_transcripts.fa, ../dbases/primary_transcripts.gtf
    output: ../dbases/salmon/hash.bin
    log: ../dbases/logs/generate_primary_transcript_salmon_index.txt
    reason: Missing output files: ../dbases/salmon/hash.bin; Input files updated by another job: ../dbases/primary_transcripts.fa, ../dbases/primary_transcripts.gtf
    resources: all_threads=12

      salmon index         -p 1         -t ../dbases/primary_transcripts.fa         --type "quasi"         -i ../dbases/salmon

rule salmon_primary:
    input: ../data/raw_data/test.fastq.gz, ../dbases/salmon/hash.bin, ../dbases/primary_transcripts.gtf
    output: ../data/salmon/test/quant.sf
    log: ../data/salmon/logs/test
    reason: Missing output files: ../data/salmon/test/quant.sf; Input files updated by another job: ../dbases/salmon/hash.bin, ../dbases/primary_transcripts.gtf
    resources: all_threads=12

      salmon quant         -i ../dbases/salmon         -r ../data/raw_data/test.fastq.gz         -o ../data/salmon/test         -p 1         --numBootstraps 50         --libType 'SR'         -g ../dbases/primary_transcripts.gtf

localrule all:
    input: ../data/salmon/test/quant.sf
    reason: Input files updated by another job: ../data/salmon/test/quant.sf
Job counts:
    count   jobs
    1   all
    1   primary_transcript_gtf
    1   salmon_idx
    1   salmon_primary
    4
```

to actually run all of these commands we submit 1 job to tesla

```bash
bsub < snakecharmer.sh
```

this script will run the snakemake pipeline. By default the script will only use up to 24 cores at a time. 
You can modify this limit if necessary by editing the following in `snakecharmer.sh`

Change jobs # and the all_threads argument to modify the number of cores or jobs run at one time. 
```bash
snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 24 \
    --resources all_threads=24 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config.yaml
```

