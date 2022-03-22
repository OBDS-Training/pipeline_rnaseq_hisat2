# OBDS-Training / pipeline_rnaseq_hisat2

A minimal pipeline demonstrating the processing paired-end RNA-sequencing using [cgatcore][link-cgatcore] and [HISAT2][link-hisat2].

## Usage

### Setup - Files

Clone this repository.

```
git clone git@github.com:OBDS-Training/pipeline_rnaseq_hisat2.git
```

Create sub-directories for the input files.

```
mkdir data
mkdir data/fastq
mkdir data/hisat2_index
```

Download an example set of compressed FASTQ files into the `data/fastq/` directory.

```
wget \
  https://github.com/sims-lab/simulated_ngs_datasets/raw/files/human.chr22.genes2/outputs/simulated_reads/sample_01_R1.fastq.gz \
  https://github.com/sims-lab/simulated_ngs_datasets/raw/files/human.chr22.genes2/outputs/simulated_reads/sample_01_R2.fastq.gz \
  https://github.com/sims-lab/simulated_ngs_datasets/raw/files/human.chr22.genes2/outputs/simulated_reads/sample_02_R1.fastq.gz \
  https://github.com/sims-lab/simulated_ngs_datasets/raw/files/human.chr22.genes2/outputs/simulated_reads/sample_02_R2.fastq.gz \
  -P data/fastq \
  --no-verbose
```

Download an example genome FASTA file into the `data/` directory, and decompress it.

```
wget \
  http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz \
  -P data \
  --no-verbose
gunzip data/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
```

Download an example GTF file into the `data/` directory.

```
wget \
  https://raw.githubusercontent.com/sims-lab/simulated_ngs_datasets/files/human.chr22.genes2/outputs/chr22.genes2.gtf \
  -P data \
  --no-verbose
```

### Setup - Environment

Create the Conda environment expected by the pipeline, using the `conda.yml`.

```
mamba env create -f conda.yml
```

Activate the Conda environment.

```
conda activate pipeline_rnaseq_hisat2
```

## Execution

Launch the pipeline.

```
python pipeline.py make full -v 5
```

The example above request the execution of entire pipeline (`make full`),
with maximal verbosity (`-v 5`).

[link-cgatcore]: https://github.com/cgat-developers/cgat-core
[link-hisat2]: http://www.ccb.jhu.edu/software/hisat/index.shtml
