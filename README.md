[![CI](https://github.com/OBDS-Training/pipeline_rnaseq_hisat2/actions/workflows/build.yml/badge.svg)](https://github.com/OBDS-Training/pipeline_rnaseq_hisat2/actions/workflows/build.yml)

# pipeline_rnaseq_hisat2

A minimal pipeline demonstrating the processing paired-end RNA-sequencing using [cgatcore][link-cgatcore] and [HISAT2][link-hisat2].

## Usage

1. Create a new repository from this one, using the `Use as template` button on [GitHub](https://github.com/OBDS-Training/pipeline_rnaseq_hisat2).
    + That way, your new repository starts its own commit history, where you can record your own changes!
    + Only fork this repository if you wish to contribute updates to the template pipeline itself.
2. Clone the new repository to the computer where you wish to run the pipeline.
    + A clone is a working directory for one run of the pipeline on one set of FASTQ files.
    + To run the pipeline on another set of FASTQ file, go back to step 1, and create another
      repository (with a different name) from the template.
3. Create a Conda environment named `pipeline_rnaseq_hisat2` using the file `envs/pipeline.yml`. 
    + You only need to do this once, no matter how many times you run the pipeline and how many
      copies of the pipeline you have cloned.
    + A Conda environment can be shared by any number of pipelines.
    + In doubt - if there is a chance that the environment was altered - remove the existing
      environment and create it again from this file.
4. Edit the set of input files in `config/fastq_links.tsv`:
    + In the first column, declare the location of input files (preferably as absolute paths,
      or relative to the working directory of the clone).
    + In the second column, declare a unique name that will be used to create a symbolic link
      to the corresponding input file.
    + Input files should be stored outside the working directory of the clone.
    + Symbolic links will be automatically created in a sub-directory `data/fastq`,
      in the working directory of the clone.
5. Edit the configuration of the pipeline as needed, in the file `config/pipeline.yml`.
    + Commit your changes to the configuration for version control and traceability.
6. Run the pipeline!
    + On a High-Performance Computing (HPC) cluster, `python pipeline.py make full -v 5`, to use the Distributed Resource Management Application API (DRMAA).
    + On a local machine `python pipeline.py make full -v 5 --local`.

[link-cgatcore]: https://github.com/cgat-developers/cgat-core
[link-hisat2]: http://www.ccb.jhu.edu/software/hisat/index.shtml
