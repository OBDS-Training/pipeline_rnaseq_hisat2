[![CI](https://github.com/sims-lab/pipeline_rnaseq_hisat2/actions/workflows/build.yml/badge.svg)](https://github.com/sims-lab/pipeline_rnaseq_hisat2/actions/workflows/build.yml)

# pipeline_rnaseq_hisat2

Pipeline for processing paired-end RNA-sequencing using [cgatcore][link-cgatcore] and [HISAT2](http://www.ccb.jhu.edu/software/hisat/index.shtml).

## Usage

1. Create a new repository from this one, using the `Use as template` button on [GitHub](https://github.com/sims-lab/pipeline_rnaseq_hisat2).
    + That way, your new repository starts its own commit history, where you can record your own changes!
    + Only fork this repository if you wish to contribute updates to the template pipeline itself.
2. Clone the new repository to the computer where you wish to run the pipeline.
    + The clone is the working directory for one run of the pipeline on one set of FASTQ files.
    + To run the pipeline on another set of FASTQ file, go back to step 1, and create another repository from the template.
3. Create a Conda environment named `pipeline_rnaseq_hisat2` using the file `envs/pipeline.yml`. 
    + You only need to do this once, no matter how many times you run the pipeline and how many copies of the pipeline you have cloned.
    + In doubt, remove the existing environment and create it again from this file.
4. Create symbolic links to your input FASTQ files, in the subdirectory `data/`.
    + Do not copy the files themselves, or make sure you don't commit them to Git (e.g. use `.gitignore`).
5. Edit the configuration of the pipeline as needed, in the file `config.yml`.
    + Commit your changes to the configuration for version control and traceability.
6. Run the pipeline!
    + On a High-Performance Computing (HPC) cluster, `python pipeline.py make full -v 5`, to use the Distributed Resource Management Application API (DRMAA).
    + On a local machine `python pipeline.py make full -v 5 --local`.

[link-cgatcore]: https://github.com/cgat-developers/cgat-core
